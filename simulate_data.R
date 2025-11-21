###########################################################################
# Generate data with the following variables:
# Z: treatment/exposure (binary)
# L: mediator-outcome confounder (binary)
# clrM: clr transformation of compositional mediators
# Y: outcome (binary)

# Main function: simulate_data
# Input: 
# seed: seed for generating random values
# n_subj: sample size
# X: baseline covariates 
# alpha_0, alpha_x: parameters for the model to generate Z
# gamma_0, gamma_z, gamma_x: parameters for the model to generate L
# eta_0, theta_0, theta_z, theta_x, theta_l, exp_order, eta_lower, eta_upper: parameters for the model to generate clrM
# beta_0, beta_z, beta_l, beta_x, beta_M: parameters for the model to generate Y

# Output:
# fulldf: a generated data frame
###########################################################################

#define subfunctions which will be used by the main function
expit <- function(x){1 / (1 + exp(-x))}
comp <- function(x){exp(x)/sum(exp(x))}

#the main function for generating the data
simulate_data <- function(seed, n_subj, X, alpha_0, alpha_x, gamma_0, gamma_z, gamma_x, eta_0, theta_0, theta_z, theta_x, theta_l, exp_order, eta_lower, eta_upper, 
                           beta_0, beta_z, beta_l, beta_x, beta_M){
  
  set.seed(seed)
  
  #generate Z
  p_z <- expit(alpha_0 + X %*% alpha_x)
  Z <-  rbinom(n_subj, 1, p_z)

  #generate L
  p_l <- expit(gamma_0 + gamma_z*Z + X %*% gamma_x)
  L <- rbinom(n_subj, 1, p_l)

  #generate clrM
  shape.val <- eta_0
  mean.pred <- theta_0 + theta_z*Z + X %*% theta_x + L*theta_l
  scale.val <- (exp(mean.pred))^exp_order / eta_0
  
  constant <- rgamma(n = n_subj, shape = shape.val, scale = scale.val)
  n.constant <- rdunif(n_subj, b = floor(eta_upper*(n_M-1)), a = floor(eta_lower*(n_M-1)))
 
  clrM = matrix(0, n_subj, n_M)
  for (i in 1:n_subj){
    set.seed((seed-1)*n_subj+i)
    clrM[i,1:(n_M-1)] <- - constant[i]
    randomvalue <- c(runif(n_M-1, 0, 1), 2)
    clrM[i,(rank(randomvalue)<(n_M-n.constant[i]))] = runif((n_M-n.constant[i]-1), - constant[i]*(1-n.constant[i])/(n_M-1-n.constant[i]), (1+n.constant[i])*constant[i]/(n_M-1-n.constant[i]))
  }  
  
  clrM[,n_M] <- - apply(clrM[,1:(n_M-1)],1,sum)
  colnames(clrM) <- paste0("clrM", seq(1:n_M))

  #generate Y  
  linear.Y <-  beta_0 + beta_z*Z + beta_l*L + X %*% beta_x + clrM %*% beta_M 
  p_y <-  expit(linear.Y)
  Y <-  rbinom(n = n_subj, 1, p_y) 
  print(paste("mean(Y) =", mean(Y)))
  
  
  fulldf <- as.data.frame(cbind(Y, Z, L, X, clrM))
  
  return(fulldf)
  
}

