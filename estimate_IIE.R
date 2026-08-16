###############################################################################
# Estimate the interventional indirect effect (IIE) using inverse probability
# weighting (IPW).
#
# Main function: IIE.IPW.est
#
# Arguments:
# temp_df:  data frame containing Y, Z, L, X, and the mediators
# beta_hat: estimated coefficients from the outcome model
#
# Value:
# Estimated IIE, calculated as eta_3 - eta_1.
###############################################################################


# Define helper functions used by the main function
expit <- function(x) {
  1 / (1 + exp(-x))
}


# Estimate the IIE using IPW
IIE.IPW.est <- function(temp_df, beta_hat) {
  # Extract variables and dimensions
  Y <- temp_df[, 1]
  Z <- temp_df[, 2]
  L <- temp_df[, 3]
  
  X <- as.matrix(temp_df[, 4:5])
  X1 <- X[, 1]
  X2 <- X[, 2]
  
  beta_M <- beta_hat[-(1:5)]
  n_M <- length(beta_M)
  
  M <- as.matrix(temp_df[, 6:(5 + n_M)])
  M_select <- as.matrix(M[, beta_M != 0, drop = FALSE])
  beta_M_select <- beta_M[beta_M != 0]
  
  # Estimate the model for L conditional on Z and X
  mle1 <- glm(
    L ~ Z + X,
    family = binomial(link = "logit")
  )
  
  gamma1_0 <- mle1$coefficients[1]
  gamma1_z <- mle1$coefficients[2]
  gamma1_x <- mle1$coefficients[3:4]
  
  linear_pred_l_z1 <- gamma1_0 +
    gamma1_z +
    X %*% gamma1_x
  
  p_l_z1 <- ifelse(
    L == 1,
    expit(linear_pred_l_z1),
    1 - expit(linear_pred_l_z1)
  )
  
  # Estimate the model for L conditional on Z, X, and selected mediators
  set.seed(1)
  
  esCtrl <- list(
    n.hidden = 5,
    activate = "relu",
    l1.reg = 0.00001,
    early.stop.det = 1000,
    n.batch = 20,
    n.epoch = 100,
    learning.rate.adaptive = "adam",
    plot = FALSE
  )
  
  Z0 <- Z
  xx0 <- as.matrix(cbind(Z0, X1, X2, M_select))
  
  dnn_obj <- importDnnet(
    x = xx0,
    y = factor(L, levels = c(0, 1))
  )
  
  fit2 <- ensemble_dnnet(
    dnn_obj,
    20,
    esCtrl,
    best.opti = TRUE,
    verbose = FALSE
  )
  
  # Estimate P(L | Z = 0, M, X)
  n_subj <- length(Y)
  Z0 <- rep(0, n_subj)
  xx1 <- as.matrix(cbind(Z0, X1, X2, M_select))
  pred1 <- predict(fit2, xx1)[, "1"]
  
  p_l_z0m <- ifelse(
    L == 1,
    pred1,
    1 - pred1
  )
  p_l_z0m <- pmax(p_l_z0m, 1e-5)
  
  # Estimate P(L | Z = 1, M, X)
  Z0 <- rep(1, n_subj)
  xx2 <- as.matrix(cbind(Z0, X1, X2, M_select))
  pred2 <- predict(fit2, xx2)[, "1"]
  
  p_l_z1m <- ifelse(
    L == 1,
    pred2,
    1 - pred2
  )
  p_l_z1m <- pmax(p_l_z1m, 1e-5)
  
  # Estimate the IIE as eta_3 - eta_1
  #
  # eta_3 =
  # E[Z Y P(L | Z = 1, X) /
  #   {P(Z = 1 | X) P(L | Z = 1, M^(1), X)}]
  #
  # eta_1 =
  # E[(1 - Z) Y E(Y | Z = 1, L, M^(1), X) P(L | Z = 1, X) /
  #   {P(Z = 0 | X) E(Y | Z = 0, L, M^(1), X)
  #    P(L | Z = 0, M^(1), X)}]
  
  # Estimate eta_3
  exposure_model <- glm(
    Z ~ X1 + X2,
    family = binomial(link = "logit")
  )
  
  p_z_hat <- predict(
    exposure_model,
    type = "response"
  )
  
  w_1 <- p_l_z1 / (p_z_hat * p_l_z1m)
  eta_3 <- mean(Z * Y * w_1)
  
  # Extract coefficients from the outcome model
  beta_0 <- beta_hat[1]
  beta_z <- beta_hat[2]
  beta_l <- beta_hat[3]
  beta_x1 <- beta_hat[4]
  beta_x2 <- beta_hat[5]
  
  # Estimate conditional outcome probabilities
  E_Y_z1_l <- expit(
    beta_0 +
      beta_z +
      beta_l * L +
      beta_x1 * X1 +
      beta_x2 * X2 +
      M_select %*% beta_M_select
  )
  
  E_Y_z0_l <- expit(
    beta_0 +
      beta_l * L +
      beta_x1 * X1 +
      beta_x2 * X2 +
      M_select %*% beta_M_select
  )
  
  # Estimate eta_1
  w_2 <- E_Y_z1_l * p_l_z1 /
    ((1 - p_z_hat) * E_Y_z0_l * p_l_z0m)
  
  eta_1 <- mean((1 - Z) * Y * w_2)
  
  # Estimate the IIE
  IIE.est <- eta_3 - eta_1
  
  return(IIE.est)
}

