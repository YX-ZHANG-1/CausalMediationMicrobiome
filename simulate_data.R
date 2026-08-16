###############################################################################
# Simulate a data set containing the following variables:
#
# X:    baseline covariates
# Z:    binary  exposure
# L:    binary exposure-induced mediator-outcome confounder
# clrM: centered log-ratio (clr) transformed mediators
# Y:    binary outcome
#
# Main function: simulate_data
#
# Arguments:
# seed:      random-number seed
# n_subj:    number of subjects
# n_M:       number of mediator components
# X:         matrix or data frame of baseline covariates
#
# alpha_0, alpha_x:
#   Parameters of the model used to generate Z.
#
# gamma_0, gamma_z, gamma_x:
#   Parameters of the model used to generate L.
#
# eta_0, theta_0, theta_z, theta_x, theta_l, eta_lower, eta_upper:
#   Parameters of the model used to generate clrM.
#
# beta_0, beta_z, beta_l, beta_x, beta_M:
#   Parameters of the model used to generate Y.
#
# Value:
# A data frame containing Y, Z, L, X, and clrM.
###############################################################################

# Define helper functions used by the main function
expit <- function(x) {
  1 / (1 + exp(-x))
}


# Generate the simulated data
simulate_data <- function(
    seed,
    n_subj,
    n_M,
    X,
    alpha_0,
    alpha_x,
    gamma_0,
    gamma_z,
    gamma_x,
    eta_0,
    theta_0,
    theta_z,
    theta_x,
    theta_l,
    eta_lower,
    eta_upper,
    beta_0,
    beta_z,
    beta_l,
    beta_x,
    beta_M
) {
  set.seed(seed)
  
  # Generate Z
  p_z <- expit(alpha_0 + X %*% alpha_x)
  Z <- rbinom(n_subj, 1, p_z)
  
  # Generate L
  p_l <- expit(gamma_0 + gamma_z * Z + X %*% gamma_x)
  L <- rbinom(n_subj, 1, p_l)
  
  # Generate clrM
  shape.val <- eta_0
  mean.pred <- theta_0 +
    theta_z * Z +
    X %*% theta_x +
    L * theta_l
  scale.val <- exp(mean.pred) / eta_0
  
  constant <- -rgamma(
    n = n_subj,
    shape = shape.val,
    scale = scale.val
  )
  
  n.constant <- sample(
    floor(eta_lower * (n_M - 1)):floor(eta_upper * (n_M - 1)),
    size = n_subj,
    replace = TRUE
  )
  
  clrM <- matrix(0, nrow = n_subj, ncol = n_M)
  
  for (i in seq_len(n_subj)) {
    set.seed((seed - 1) * n_subj + i)
    
    clrM[i, 1:(n_M - 1)] <- constant[i]
    
    randomvalue <- runif(n_M - 1)
    idM <- which(
      rank(randomvalue) <= (n_M - 1 - n.constant[i])
    )
    
    clrM[i, idM] <- runif(
      length(idM),
      constant[i] * (1 - n.constant[i]) /
        (n_M - 1 - n.constant[i]),
      -(1 + n.constant[i]) * constant[i] /
        (n_M - 1 - n.constant[i])
    )
  }
  
  clrM[, n_M] <- -apply(
    clrM[, 1:(n_M - 1), drop = FALSE],
    1,
    sum
  )
  colnames(clrM) <- paste0("clrM", seq_len(n_M))
  
  # Generate Y
  linear.Y <- beta_0 +
    beta_z * Z +
    beta_l * L +
    X %*% beta_x +
    clrM %*% beta_M
  
  p_y <- expit(linear.Y)
  Y <- rbinom(n = n_subj, size = 1, prob = p_y)
  
  print(paste("mean(Y) =", mean(Y)))
  
  fulldf <- as.data.frame(cbind(Y, Z, L, X, clrM))
  
  return(fulldf)
}
