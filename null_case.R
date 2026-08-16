###############################################################################
# Monte Carlo simulation for estimating the interventional indirect effect
# under a null hypothesis with a zero mediation effect.
###############################################################################


# Clear the current R environment
rm(list = ls())

# Record the starting time
ptm <- proc.time()

# Set the working directory to the directory containing this script.
# Alternatively, run this script from the project root and remove setwd().
# setwd("/path/to/project")

# Load functions
source("simulate_data.R")
source("estimate_IIE.R")
source("Outcome_model_fit.R")
source("bootstrap.R")

# Load baseline covariates
realdata <- as.matrix(
  read.csv(
    "a pseudo-dataset for baseline covariates X.csv",
    header = TRUE,
    sep = ","
  )
)

realX <- data.frame(realdata)

# The file above is a pseudo-dataset used in place of the real baseline
# covariate data.


###############################################################################
# Simulation settings
###############################################################################

# Simulation case
case <- "null"

# Number of Monte Carlo replications
nrep <- 500

# Number of bootstrap replications
nbts <- 400

# Sample size for each Monte Carlo data set
n_subj <- 70

# Number of CPUs used for parallel computing
ncpus <- 40

# Parameters for the exposure model
alpha_0 <- 0.2
alpha_x <- c(-0.5, 0.5)

# Parameters for the exposure-induced mediator-outcome confounder model
gamma_0 <- 0.5
gamma_z <- 0
gamma_x <- c(0.5, -0.5)

# Parameters for the centered log-ratio (clr) transformed mediator model
n_M <- 134
eta_0 <- 4
theta_0 <- -1
theta_z <- 0
theta_x <- c(-0.8, -0.2)
theta_l <- 0

eta_lower <- 0.4
eta_upper <- 0.8

# Parameters for the outcome model
beta_0 <- 3
beta_z <- -1.2
beta_l <- -6
beta_x <- as.matrix(c(-1.2, -1.2))
beta_M <- rep(0, n_M)

beta.true <- c(
  beta_0,
  beta_z,
  beta_l,
  beta_x,
  beta_M
)

# Variable-selection and tuning methods
penalty <- "SCAD"
tune <- "aic"

# Maximum number of estimation attempts
iter1 <- 50

# Maximum acceptable absolute value of an IIE estimate
threshold <- 1

# Save the simulation settings
true.parameters <- list(
  nrep = nrep,
  nbts = nbts,
  n_subj = n_subj,
  alpha_0 = alpha_0,
  alpha_x = alpha_x,
  gamma_0 = gamma_0,
  gamma_z = gamma_z,
  gamma_x = gamma_x,
  n_M = n_M,
  eta_0 = eta_0,
  theta_0 = theta_0,
  theta_z = theta_z,
  theta_x = theta_x,
  theta_l = theta_l,
  eta_lower = eta_lower,
  eta_upper = eta_upper,
  beta.true = beta.true,
  penalty = penalty,
  tune = tune,
  iter1 = iter1,
  threshold = threshold
)

# True value of IIE 
IIE.true <- 0


###############################################################################
# Run one Monte Carlo replication
###############################################################################

wrapper <- function(s) {
  condition <- TRUE
  counter <- 0
  success_b <- FALSE
  
  IIE_attempt <- list(
    success = FALSE,
    val = NULL
  )
  
  while ((!success_b || condition) && counter < iter1) {
    # Generate a seed for the current estimation attempt
    seed <- (s - 1) * iter1 + counter
    set.seed(seed)
    
    # Generate X by resampling the baseline covariate data
    X0 <- sample_n(
      realX,
      size = n_subj,
      replace = TRUE
    )
    
    X <- matrix(
      as.numeric(as.matrix(X0)),
      ncol = 2
    )
    
    X[, 1] <- X[, 1] / 100
    
    # Simulate the analysis data
    simu_df <- simulate_data(
      seed = seed,
      n_subj = n_subj,
      n_M = n_M,
      X = X,
      alpha_0 = alpha_0,
      alpha_x = alpha_x,
      gamma_0 = gamma_0,
      gamma_z = gamma_z,
      gamma_x = gamma_x,
      eta_0 = eta_0,
      theta_0 = theta_0,
      theta_z = theta_z,
      theta_x = theta_x,
      theta_l = theta_l,
      eta_lower = eta_lower,
      eta_upper = eta_upper,
      beta_0 = beta_0,
      beta_z = beta_z,
      beta_l = beta_l,
      beta_x = beta_x,
      beta_M = beta_M
    )
    
    # Count the number of estimation attempts
    counter <- counter + 1
    
    # Prepare the outcome and predictors
    x <- as.matrix(simu_df[, -1])
    y <- simu_df[, 1]
    
    # Do not penalize Z, L, X1, or X2
    penalty.factor <- c(
      rep(0, 4),
      rep(1, n_M)
    )
    
    # Estimate the outcome model
    model_attempt <- tryCatch(
      {
        model <- tune.fit(
          x,
          y,
          family = "binomial",
          penalty = penalty,
          type.measure = "deviance",
          tune = tune,
          nfolds = 5,
          penalty.factor = penalty.factor
        )
        
        beta_hat <- as.matrix(
          model$coef.beta
        )
        
        invalid_beta <- is.null(beta_hat) ||
          length(beta_hat) == 0 ||
          any(is.na(beta_hat)) ||
          any(!is.finite(beta_hat)) ||
          all(beta_hat == 0)
        
        if (invalid_beta) {
          stop("Invalid estimator of beta for IIE.")
        }
        
        list(
          success = TRUE,
          val = beta_hat
        )
      },
      error = function(e) {
        list(
          success = FALSE,
          val = NULL
        )
      }
    )
    
    if (!model_attempt$success) {
      next
    }
    
    beta_hat <- model_attempt$val
    
    # Estimate the IIE
    IIE_attempt <- tryCatch(
      {
        IIE_hat <- IIE.IPW.est(
          temp_df = simu_df,
          beta_hat = beta_hat
        )
        
        list(
          success = TRUE,
          val = IIE_hat
        )
      },
      error = function(e) {
        list(
          success = FALSE,
          val = NULL
        )
      }
    )
    
    if (IIE_attempt$success && !is.null(IIE_attempt$val)) {
      success_b <- TRUE
      IIE_hat <- IIE_attempt$val
      
      # Regenerate the data if the estimate is missing or exceeds the
      # prespecified threshold.
      condition <- is.na(IIE_hat) ||
        abs(IIE_hat) > threshold
    } else {
      success_b <- FALSE
      condition <- TRUE
    }
  }
  
  if (success_b && !condition) {
    IIE_hat <- IIE_attempt$val
    
    # Obtain bootstrap estimates of the IIE
    IIE_hat.b <- IPW.bootstrap(
      data = simu_df,
      nbts = nbts,
      n_M = n_M,
      penalty = penalty,
      tune = tune,
      iter1 = iter1,
      threshold = threshold
    )
  } else {
    IIE_hat <- NA_real_
    
    IIE_hat.b <- tibble(
      IIE_bs_all = rep(NA_real_, nbts)
    )
  }
  
  return(
    list(
      IIE_hat = IIE_hat,
      IIE_hat.b = IIE_hat.b
    )
  )
}


###############################################################################
# Parallel computation
###############################################################################

library(snowfall)

sfInit(
  parallel = TRUE,
  cpus = ncpus
)

sfExportAll(
  except = NULL,
  debug = FALSE
)

sfLibrary(deepTL)
sfLibrary(MASS)
sfLibrary(dplyr)
sfLibrary(tidyr)
sfLibrary(tibble)
sfLibrary(dglm)
sfLibrary(ncvreg)
sfLibrary(glmnet)



outputs <- sfLapply(
  seq_len(nrep),
  wrapper
)


###############################################################################
# Combine the Monte Carlo and bootstrap results
###############################################################################

IIE_hat_all <- rep(NA_real_, nrep)
IIE_hat_bs_all <- vector(
  mode = "list",
  length = nrep
)
IIE.bs.all <- vector(
  mode = "list",
  length = nrep
)

for (j in seq_len(nrep)) {
  IIE_hat_all[j] <- outputs[[j]]$IIE_hat
  IIE_hat_bs_all[[j]] <- outputs[[j]]$IIE_hat.b
  
  IIE.bs.all[[j]] <- IIE_hat_bs_all[[j]] %>%
   mutate(seed = j)
}


###############################################################################
# Summarize the simulation results
###############################################################################

# Monte Carlo bias
IIE_bias <- mean(
  IIE_hat_all - IIE.true,
  na.rm = TRUE
)

# Monte Carlo standard deviation
IIE_sd <- sd(
  IIE_hat_all,
  na.rm = TRUE
)

# Root mean squared error
IIE_rmse <- sqrt(
  mean(
    (IIE_hat_all - IIE.true)^2,
    na.rm = TRUE
  )
)

# Original IIE estimates
IIE_hat.ori <- tibble(
  IIE_hat = IIE_hat_all,
  seed = seq_len(nrep)
)

# Construct bootstrap confidence intervals
IIE_tbl <- do.call("rbind", IIE.bs.all) %>%
  mutate(
    IIE_bs_all = as.numeric(IIE_bs_all),
    seed = as.integer(seed)
  ) %>%
  left_join(
    IIE_hat.ori,
    by = "seed"
  ) %>%
  group_by(
    seed,
    IIE_hat
  ) %>%
  drop_na() %>%
  summarise(
    IIE_b_mean = mean(
      IIE_bs_all,
      na.rm = TRUE
    ),
    sd = sd(
      IIE_bs_all,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Normal.lower = IIE_hat +
      qnorm(0.025) * sd,
    Normal.upper = IIE_hat +
      qnorm(0.975) * sd,
    Normal.cover = (
      Normal.lower < 0 &
        Normal.upper > 0
    ),
    case = case,
    n = n_subj
  )

# Calculate the type I error
type_one_error <- IIE_tbl %>%
  group_by(case, n) %>%
  summarise(
    Normal.CI = 1 - mean(
      Normal.cover,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


###############################################################################
# Finish parallel processing and save the results
###############################################################################

elapsed_time <- proc.time()[3] - ptm[3]
print(elapsed_time)

sfStop()

output_file <- paste0(
  "./results_penalty_",
  penalty,
  "_tune_",
  tune,
  "_case",
  case,
  "_n",
  n_subj,
  "_nrep",
  nrep,
  ".RData"
)

save(
  true.parameters,
  IIE.true,
  IIE_hat_all,
  IIE_hat_bs_all,
  IIE_bias,
  IIE_sd,
  IIE_rmse,
  type_one_error,
  elapsed_time,
  file = output_file
)
