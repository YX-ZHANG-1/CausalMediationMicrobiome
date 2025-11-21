###########################################################################
# Estimate the true value of IIE based on a generated large sample

# Main function: IIE.IPW.true
# Input: 
# n_subj: sample size
# seed: seed for generating random values
# X: baseline covariates
# alpha_0, alpha_x: parameters for the model to generate Z
# gamma_0, gamma_z, gamma_x: parameters for the model to generate L
# eta_0, theta_0, theta_z, theta_x, theta_l, exp_order, eta_lower, eta_upper: parameters for the model to generate clrM
# beta_0, beta_z, beta_l, beta_x, beta_M: parameters for the model to generate Y

# Output:
# IIE: estimated value of true IIE 
###########################################################################

#define a subfunction which will be used by the other function
expit <- function(x){1 / (1 + exp(-x))}


#the main function for estimating the true value of IIE
IIE.IPW.true <- function(n_true, seed, X, alpha_0, alpha_x, gamma_0, gamma_z, gamma_x, eta_0, theta_0, theta_z, theta_x, theta_l, exp_order, eta_lower, eta_upper, 
                         beta_0, beta_z, beta_l, beta_x, beta_M){
  
  #generate data
  temp_df <- simulate_data(seed, n_true, X, alpha_0, alpha_x, gamma_0, gamma_z, gamma_x, eta_0, theta_0, theta_z, theta_x, theta_l, exp_order, eta_lower, eta_upper, 
                           beta_0, beta_z, beta_l, beta_x, beta_M)
  
  Y <- temp_df[,1]
  Z <- temp_df[,2]
  L <- temp_df[,3]
  X <- as.matrix(temp_df[,4:5])
  X1 <- X[,1]
  X2 <- X[,2]
  n_M <- length(beta_M)
  clrM <- as.matrix(temp_df[,6:(5+n_M)])
  clrM.min <- apply(clrM, 1, min)


  #calculate values of P(l | z = 1, x) and P(l | z = 0, x)
  p_l_z1 <-  ifelse(L == 1, expit(gamma_0 + gamma_z + X %*% gamma_x), 1 - expit(gamma_0 + gamma_z + X %*% gamma_x))
  p_l_z0 <-  ifelse(L == 1, expit(gamma_0 + X %*% gamma_x), 1 - expit(gamma_0 + X %*% gamma_x))
  
    
  #calculate P(clrm | z, l, x)
  b = floor(eta_upper * (n_M - 1))
  a = floor(eta_lower * (n_M - 1))
  
  
  f_clrm1m2_zl <- function(data){
    m1 <- data[1]
    m2 <- data[2]
    z <- data[3]
    l <- data[4]
    x <- data[5:6]
    m.min <- data[7]
    
    scale <- (exp(theta_0 + theta_z * z + theta_l * l + x %*% theta_x))^exp_order / eta_0
    
    if (m1 == m2){
      fm1m2_zlx <- (((b - a + 1)^2 - 1) / 12 + ((a + b) / 2)^2 - (a + b) / 2) / ((n_M - 1) * (n_M - 2)) * dgamma(-m1, shape = eta_0, scale = scale)
      
    }else if(m1 > m2){
      if(m2 == m.min){
        
        part10 <- function(s){
          1 / (-(1 + s) * m2 / (n_M - 1 - s) - m2*(1 - s) / (n_M - 1 - s)) * (s / (n_M - 1)) * (1 - s / (n_M - 1)) / (b - a + 1) * dgamma(-m2, shape = eta_0, scale = scale)
        } 
        
        fm1m2_zlx <- sum(sapply(a:b, part10))
        
      }else{
        
        part20 <- function(c){
          1 / (-(1 + s) * c / (n_M - 1 - s) - m2) / (-(1 + s) * c / (n_M - 1 - s) - c * (1 - s) / (n_M - 1 - s)) * ((n_M - 1 - s) * (n_M - 2 - s) / ((n_M - 1) * (n_M - 2))) / 2 / (b - a + 1) * dgamma(-c, shape = eta_0, scale = scale)
        }
        
        part2 <- 0
        
        for (s in a:b){
          upperbound <- min(0, -m1 * (n_M - 1 - s) / (1 + s))
          lowerbound <- m2 * (n_M - 1 - s) / (1 - s)
          part21 <- cubintegrate(f = part20, lower = lowerbound, upper = upperbound, method = "hcubature")$integral
          part2 <- part2 + part21
        }
        
        fm1m2_zlx <- part2
      }
      
    }else{
      
      if(m1 == m.min){
        
        part12 <- function(s){
          1 / (-(1 + s) * m1 / (n_M - 1 - s) - m1 * (1 - s) / (n_M - 1 - s)) * (s / (n_M - 1)) * (1 - s / (n_M - 1)) / (b - a + 1) * dgamma(-m1, shape = eta_0, scale = scale)
        } 
        
        fm1m2_zlx <- sum(sapply(a:b, part12))
        
      }else{
        
        part22 <- function(c){
          1/(-(1 + s) * c/(n_M - 1 - s) - m1)/(-(1 + s) * c / (n_M - 1 - s) - c* (1 - s) / (n_M - 1 - s)) * ((n_M - 1 - s) * (n_M - 2 - s) / ((n_M - 1) * (n_M - 2))) / 2 / (b - a + 1) * dgamma(-c, shape = eta_0, scale = scale)
        }
        
        part2 <- 0
        
        for (s in a:b){
          upperbound <- min(0, -m2 * (n_M - 1 - s) / (1 + s))
          lowerbound <- m1 * (n_M - 1 - s) / (1 - s)
          part21 <- cubintegrate(f = part22, lower = lowerbound, upper = upperbound, method = "hcubature")$integral
          part2 <- part2 + part21
        }
        
        fm1m2_zlx <- part2
      }
    }
    return(fm1m2_zlx) 
  }
  
  ##calculate P(clrm | z = 1, l, x)
  data1 <- as.matrix(cbind(clrM[,1:2], rep(1,n_true), L, X, clrM.min))
  f_clrm_z1l <- apply(data1,1,f_clrm1m2_zl)
  
  ##calculate P(clrm | z = 1, l = 0, x)
  data2 <- as.matrix(cbind(clrM[,1:2], rep(1,n_true), rep(0,n_true), X, clrM.min))
  f_clrm_z1l0 <- apply(data2, 1, f_clrm1m2_zl)
  
  ##calculate P(clrm | z = 1, l = 1, x)
  data3 <- as.matrix(cbind(clrM[,1:2], rep(1,n_true), rep(1,n_true), X, clrM.min))
  f_clrm_z1l1 <- apply(data3, 1, f_clrm1m2_zl)
  
  ##calculate P(clrm | z = 0, l, x)
  data4 <- as.matrix(cbind(clrM[,1:2], rep(0,n_true), L, X, clrM.min))
  f_clrm_z0l <- apply(data4, 1, f_clrm1m2_zl)
  
  ##calculate P(clrm | z = 0, l = 0, x)
  data5 <- as.matrix(cbind(clrM[,1:2], rep(0,n_true), rep(0,n_true), X, clrM.min))
  f_clrm_z0l0 <- apply(data5, 1, f_clrm1m2_zl)
  
  ##calculate P(clrm | z = 0, l = 1, x)
  data6 <- as.matrix(cbind(clrM[,1:2], rep(0,n_true), rep(1,n_true), X, clrM.min))
  f_clrm_z0l1 <- apply(data6, 1, f_clrm1m2_zl)
  
  ##calculate P(l | z = 1, clrm, x)
  p_l0_z1 <- 1 - expit(gamma_0 + gamma_z + X %*% gamma_x)
  p_l1_z1 <- 1 - p_l0_z1
  
  p_l_z1clrm <- f_clrm_z1l * p_l_z1 / (f_clrm_z1l0 * p_l0_z1 + f_clrm_z1l1 * p_l1_z1)
  
  ##calculate P(l | z = 0, clrm, x)
  p_l0_z0 <- 1 - expit(gamma_0 + X %*% gamma_x)
  p_l1_z0 <- 1 - p_l0_z0
  
  p_l_z0clrm <- f_clrm_z0l*p_l_z0/(f_clrm_z0l0*p_l0_z0 + f_clrm_z0l1*p_l1_z0)
  
  
  ###calculate IIE
  #IIE = eta_3 - eta_1
  #eta_3 = E[Z * Y * P(L | Z = 1, X) / (P(Z = 1 | x)  * P(L | Z = 1, M^{(1)}, X))] 
  #eta_1 = E[(1 - Z) * Y * E[Y | Z = 1, L, M^{(1)}, X] * P(L | Z = 1, X) / (P(Z = 0 | X) * E[Y | Z = 0, L, M^{(1)}, X] * P(L | Z = 0, M^{(1)}, X))]
  
  #calculate eta_3
  p_z <- expit(alpha_0 + X %*% alpha_x)
  w_1 <- p_l_z1 / (p_z*p_l_z1clrm)
  eta_3 <- mean(Z * Y * w_1)
  
  
  #calculate eta_1
  E_Y_z1 <- expit(beta_0 + beta_z + beta_l * L + X %*% beta_x + clrM %*% beta_M)
  E_Y_z0 <- expit(beta_0 + beta_l * L + X %*% beta_x + clrM %*% beta_M)

  w_2 <- E_Y_z1 * p_l_z1 / ((1 - p_z) * E_Y_z0 * p_l_z0clrm)
  eta_1 <- mean((1 - Z) * Y * w_2)

  #calculate IIE
  IIE <- eta_3 - eta_1
  
  return(IIE)
}


