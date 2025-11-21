###########################################################################
# Estimate the IIE using IPW approach

# Input: 
# a data frame with all variables
# estimated beta coefficients (in the Y model)

# Output:
# IIE.est: estimated value of the IIE 
###########################################################################


expit <- function(x){1 / (1 + exp(-x))}


IIE.IPW.est <- function(temp_df, beta_hat){
  
  Y <- temp_df[,1]
  n_subj <- length(Y)
  Z <- temp_df[,2]
  L <- temp_df[,3]
  X <- as.matrix(temp_df[,4:5])
  X1 <- X[,1]
  X2 <- X[,2]
  beta_M <- beta_hat[-c(1:5)]
  n_M <- length(beta_M)
  M <- as.matrix(temp_df[,6:(5 + n_M)])
  M_select <- as.matrix(M[,(beta_M != 0)])
  beta_M_select <- beta_M[beta_M != 0]
  

  #estimate the model of L given Z, X
  mle1 <- glm(L ~ Z + X, family = binomial(link = "logit"))
  gamma1_0 <- mle1$coefficients[1]
  gamma1_z <- mle1$coefficients[2]
  gamma1_x <- mle1$coefficients[3:4]
  
  p_l_z1 <- ifelse(L == 1, expit(gamma1_0 + gamma1_z + X %*% gamma1_x), 1 - expit(gamma1_0 + gamma1_z + X %*% gamma1_x))
  

  #estimate the model of L given Z, X and selected mediators
  set.seed(1)
  esCtrl <- list(n.hidden = c(30, 30, 30), activate = "relu",
                 l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
                 n.epoch = 100, learning.rate.adaptive = "adam", plot = FALSE)
  Z0 <- Z
  xx0 <- as.matrix(cbind(Z0, X1, X2, M_select))
  dnn_obj <- importDnnet(x = xx0, y = L)
  fit2 <- ensemble_dnnet(dnn_obj, 100, esCtrl, best.opti = TRUE, verbose = FALSE)
  
  Z0 <- rep(0,n_subj)
  xx1 <- as.matrix(cbind(Z0, X1, X2, M_select))
  pred1 <- predict(fit2, xx1)
  
  p_l_z0m <- ifelse(L == 1, pred1, 1-pred1)
  p_l_z0m_adjust <- ifelse(p_l_z0m == 0, p_l_z0m + 10^(-5), p_l_z0m)
  p_l_z0m <- p_l_z0m_adjust
  
  Z0 <- rep(1,n_subj)
  xx2 <- as.matrix(cbind(Z0, X1, X2, M_select))
  pred2 <- predict(fit2, xx2)
  
  p_l_z1m <- ifelse(L == 1, pred2, 1-pred2)
  p_l_z1m_adjust <- ifelse(p_l_z1m == 0, p_l_z1m + 10^(-5), p_l_z1m)
  p_l_z1m <- p_l_z1m_adjust
  
  
  ###estimate IIE
  #IIE = eta_3 - eta_1
  #eta_3 = E[Z * Y * P(L | Z = 1, X) / (P(Z = 1 | x)  * P(L | Z = 1, M^{(1)}, X))] 
  #eta_1 = E[(1 - Z) * Y * E[Y | Z = 1, L, M^{(1)}, X] * P(L | Z = 1, X) / (P(Z = 0 | X) * E[Y | Z = 0, L, M^{(1)}, X] * P(L | Z = 0, M^{(1)}, X))]
  
  #estimate eta_3
  p_z_hat <- predict(glm(Z ~ X1 + X2, family = binomial(link = "logit")), type = "response")
  w_1 <- p_l_z1 / (p_z_hat * p_l_z1m)
  eta_3 <- mean(Z * Y * w_1)
  
  
  #estimate eta_1
  beta_0 <- beta_hat[1]
  beta_z <- beta_hat[2]
  beta_l <- beta_hat[3]
  beta_x1 <- beta_hat[4]
  beta_x2 <- beta_hat[5]
  
  E_Y_z1_l <- expit(beta_0 + beta_z + beta_l * L+ beta_x1 * X1 + beta_x2 * X2 + M_select %*% beta_M_select)
  E_Y_z0_l <- expit(beta_0 + beta_l * L + beta_x1 * X1 + beta_x2 * X2 + M_select %*% beta_M_select)

  w_2 <- E_Y_z1_l * p_l_z1 / ((1 - p_z_hat) * E_Y_z0_l * p_l_z0m)
  
  eta_1 <- mean((1 - Z) * Y * w_2)
 
   
  #estimate IIE
  IIE.est <- eta_3 - eta_1
  
  
  return(IIE.est)
}


