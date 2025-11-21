rm(list=ls())

#set the working directory based on your own address of code
setwd( )

source("simulate_data.R")
source("estimate_IIE.R")
source("Outcome_model_fit.R")
source("bootstrap.R")


#read X in real data into R
realdata <- as.matrix(read.csv("a pseudo-dataset for baseline covariates X.csv", header=TRUE, sep=","))
realX <- data.frame(realdata)
####note: "a pseudo-dataset for baseline covariates X.csv" is a pseudo-dataset in place of the real dataset for baseline covariates X.


#indicate it is the code for null case
case <- "null"

###set the parameters
#set replication times
nrep <- 500

#set the number of bootstrap replications
nbts <- 400

#set the sample size
n_subj <- 70 

#set the number of CPUs for parallel computing, the number of CPUs can be chosen based on the characteristics of cluster
ncpus <- 20

#set parameters for Z model
alpha_0 <- 0.2
alpha_x <- c(-1, 1)

#set parameters for L model
gamma_0 <- 0.5
gamma_z <- 0
gamma_x <- c(0.5, -0.5)

#set parameters for M model
n_M <- 134
eta_0 <- 4
theta_0 <- -1
theta_z <- 0
theta_x <- c(-0.8, -0.2)
theta_l <- 0
exp_order <- 1

eta_lower <- 0.4
eta_upper <- 0.8

#set parameters for Y model
beta_0 <- 3
beta_z <- -1
beta_l <- -8
beta_x <- as.matrix(c(-1, -1))
beta_M <- as.matrix(c(-8, -8, rep(0,n_M-2)))
beta.true <- c(beta_0, beta_z, beta_l, beta_x, beta_M)

#choose variable selection method
penalty <- "SCAD"
tune <- "aic"

#set the maximum iteration times for iterations
iter1 <- 50

#set the threshold value for the estimated values
threshold <- 1


true.parameters <- list(nrep=nrep, nbts=nbts, n_subj=n_subj, alpha_0=alpha_0, alpha_x=alpha_x, gamma_0=gamma_0, gamma_z=gamma_z, gamma_x=gamma_x,
                       n_M=n_M, eta_0=eta_0, theta_0=theta_0, theta_z=theta_z, theta_x=theta_x, theta_l=theta_l, exp_order=exp_order,
                       eta_lower=eta_lower, eta_upper=eta_upper, beta.true=beta.true, penalty=penalty, tune=tune)

#the true value of IIE under null case
IIE.true <- 0


ptm <- proc.time()


wrapper <- function(s){


  condition1 = TRUE
  condition2 = TRUE
  condition3 = TRUE
  counter = 0
  

  while((condition1| condition2| condition3) & (counter < iter1)){
    ###simulate data
    #generate X from the real data by resampling
    seed <- (s - 1) * iter1 + counter
    set.seed(seed)
    X0 <- sample_n(realX, size = n_subj, replace = TRUE)
    X <- matrix(as.numeric(as.matrix(X0)), ncol = 2)
    X[,1] <- X[,1]/100
    
    #generate other variables
    simu_df <- simulate_data(seed, n_subj, X, alpha_0, alpha_x, gamma_0, gamma_z, gamma_x, eta_0, theta_0, theta_z, theta_x, theta_l, exp_order, eta_lower, eta_upper, 
                             beta_0, beta_z, beta_l, beta_x, beta_M)
    
    #count the number of iterations
    counter <- counter + 1

    #estimate the Y model
    x <- as.matrix(simu_df[,-1])
    y <- simu_df[,1]
    
    optobj_M = rep(Inf, n_M)
    coef_M = matrix(0, (n_M+5), n_M)
    
    for (j in 1:n_M){
      #for each j, estimate the Y model by adding the constraint that the estimated model contains at least L,Z,X and the jth mediator
      penalty.factor0 <- c(rep(0, 4), rep(1, n_M))
      penalty.factor0[4+j] <- 0

        a <- list()
        a <- tryCatch({
          
          model <- tune.fit(x, y, family='binomial', penalty = penalty, tune=tune, nfolds = 5, type.measure = "deviance", penalty.factor=penalty.factor0)
          lambda_hat = model$lambda
          optobj_hat <- model$optobj
          beta_hat <- as.matrix(model$coef.beta)
          
          indicator <- max(abs(beta_hat - rep(0, n_M + 5)))
          
          if (indicator == 0){
            stop('the estimator of beta for IIE does not converge')
          }
          
          results <- list(lambda_hat = lambda_hat, optobj_hat = optobj_hat, beta_hat = beta_hat)
          
          list(val = results, error = NA)
        },
        
        error = function(e) { 
          list(val = NA, error = e)
        })
        
        #if there is any error in the estimation of Y model, then stop the loop
        signal <- a$error[1]
        condition1 <- (is.na(signal)==FALSE)
        
        if(condition1){
          break;
        }
      
        if (length(a$val$lambda_hat)==1){
          optobj_M[j] = a$val$optobj_hat
          coef_M[,j]= a$val$beta_hat
        }
    }
    
    
    b <- list()
    b <- tryCatch({
      
      #choose the optimal Y model which has the smallest objective function
      M.index <- which(rank(optobj_M) == (rank(optobj_M) == 1))
      beta_hat <- as.matrix(coef_M[,M.index])

      #get the estimated value of IIE
      IIE_hat <- IIE.IPW.est(temp_df = simu_df, beta_hat = beta_hat)

      
      list(val = IIE_hat, error = NA)
    },
    
    error = function(e) { 
      list(val = NA, error = e)
    })
    
    #if there is any error in the estimation of IIE, then regenerate the data
    signal2 <- b$error[1]
    condition2 <- (is.na(signal2) == FALSE)
    
    if (!condition2){
      IIE_hat <- b$val
      
      #if the estimated value of IIE is NA or its absolute value is larger than the threshold value, then regenerate the data
      condition3 <- (is.na(IIE_hat) == TRUE | abs(IIE_hat) > threshold)
    }
    
  }
  
  
  if (counter < iter1){
    IIE_hat <- b$val

    #get the results based on bootstrap samples
    IIE_hat.b <- IPW.bootstrap(data = simu_df, nbts = nbts, n_M = n_M, penalty = penalty, tune = tune, iter1 = iter1, threshold=threshold)
    
    result <- list(IIE_hat, IIE_hat.b)
  
    }else{
      
    IIE_hat <- NA
    IIE_bs_all <- rep(NA, nbts)
    IIE_hat.b <- list(IIE_bs_all = IIE_bs_all) %>% as_tibble()

    result <- list(IIE_hat, IIE_hat.b)
    
  }
  
  return(result)
}

#set to use parallel computing
library(snowfall)
sfInit(parallel = TRUE, cpus = ncpus)  

sfExportAll(except=NULL, debug=FALSE)
sfLibrary(deepTL)
sfLibrary(MASS)
sfLibrary(dplyr)
sfLibrary(tidyverse) 
sfLibrary(dglm)
sfLibrary(ncvreg)
sfLibrary(glmnet)
sfLibrary(cubature)

outputs <- sfLapply(1:nrep, wrapper)

IIE_hat_all <- rep(0, nrep)
IIE_hat_bs_all <- vector(mode = "list", length = nrep)
IIE.bs.all <- vector(mode = "list", length = nrep)


for (j in 1:nrep){
  IIE_hat_all[j] <- outputs[[j]][[1]]
  IIE_hat_bs_all[[j]] <- outputs[[j]][[2]]
  IIE.bs.all[[j]] <- IIE_hat_bs_all[[j]] %>% mutate(seed = j)
}


#analyze the results based on Monte Carlo samples
IIE_bias <- mean(IIE_hat_all - IIE.true, na.rm = TRUE)
IIE_sd <- sd(IIE_hat_all, na.rm = TRUE)
IIE_rmse <- sqrt(mean((IIE_hat_all - IIE.true)^2, na.rm = TRUE))


#analyze the results based on bootstrap samples
IIE_hat.ori <- tibble(IIE_hat = IIE_hat_all, seed = 1:nrep) 


IIE_tbl <-
  do.call("rbind", IIE.bs.all)  %>% 
  mutate(IIE_bs_all = as.numeric(IIE_bs_all),
         seed = as.integer(seed)) %>% 
  group_by(seed) %>%  
  left_join(IIE_hat.ori, by = "seed") %>% 
  group_by(seed, IIE_hat) %>% 
  drop_na() %>% 
  summarise(IIE_b_mean = mean(IIE_bs_all, na.rm = TRUE),
            sd = sd(IIE_bs_all, na.rm = TRUE)) %>% 
  mutate(Normal.lower = IIE_hat + qnorm(0.025)*sd,
         Normal.upper = IIE_hat + qnorm(0.975)*sd,
         Normal.correct = (Normal.lower < IIE.true & Normal.upper > IIE.true)) %>% 
  mutate(case = case, n = n_subj, replication = nrep) %>% 
  select(case, n, replication, everything())


#calculate the type one error
type_one_error <- IIE_tbl %>% 
  group_by(case, n, replication) %>% 
  summarise(`Normal.CI`= 1 - mean(Normal.correct, na.rm = TRUE))


#output time for running the code
proc.time()[3] - ptm[3]

sfStop()

#save the results 
save(true.parameters, IIE.true, IIE_hat_all, IIE_hat_bs_all, IIE_bias, IIE_sd, IIE_rmse, type_one_error, file=paste("./results_penalty_", penalty, "_tune_", tune, "_case", case, "_n", n_subj, "_nrep", nrep, ".RData", sep=""))

