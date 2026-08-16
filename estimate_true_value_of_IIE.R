###############################################################################
# Estimate the true interventional indirect effect (IIE) using a large
# simulated sample.
#
# Main function: IIE.IPW.true
#
# Arguments:
# seed:      random-number seed
# n_true:    sample size used to estimate the true IIE
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
# A numeric estimate of the true IIE.
###############################################################################


# Define helper functions used by the main function
expit <- function(x) {
  1 / (1 + exp(-x))
}


# Estimate the true IIE
IIE.IPW.true <- function(
    seed,
    n_true,
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
  # Generate a large simulated data set
  temp_df <- simulate_data(
    seed = seed,
    n_subj = n_true,
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
  
  # Extract variables
  Y <- temp_df[, 1]
  Z <- temp_df[, 2]
  L <- temp_df[, 3]
  X <- as.matrix(temp_df[, 4:5])
  
  clrM <- as.matrix(
    temp_df[, 6:(5 + n_M)]
  )
  clrM.min <- apply(clrM, 1, min)
  
  # Calculate P(L | Z = 1, X) and P(L | Z = 0, X)
  linear_pred_l_z1 <- gamma_0 +
    gamma_z +
    X %*% gamma_x
  
  linear_pred_l_z0 <- gamma_0 +
    X %*% gamma_x
  
  p_l_z1 <- ifelse(
    L == 1,
    expit(linear_pred_l_z1),
    1 - expit(linear_pred_l_z1)
  )
  
  p_l_z0 <- ifelse(
    L == 1,
    expit(linear_pred_l_z0),
    1 - expit(linear_pred_l_z0)
  )
  
  # Define the lower and upper numbers of constant mediator components
  a <- floor(eta_lower * (n_M - 1))
  b <- floor(eta_upper * (n_M - 1))
  
  # Calculate the joint density of the first two clrM components
  f_clrm1m2_zl <- function(data) {
    m1 <- data[1]
    m2 <- data[2]
    z <- data[3]
    l <- data[4]
    x <- data[5:6]
    m.min <- data[7]
    
    scale_val <- exp(
      theta_0 +
        theta_z * z +
        theta_l * l +
        x %*% theta_x
    ) / eta_0
    
    if (m1 == m2) {
      fm1m2_zlx <- (
        ((b - a + 1)^2 - 1) / 12 +
          ((a + b) / 2)^2 -
          (a + b) / 2
      ) /
        ((n_M - 1) * (n_M - 2)) *
        dgamma(
          -m1,
          shape = eta_0,
          scale = scale_val
        )
    } else if (m1 > m2) {
      if (m2 == m.min) {
        part10 <- function(s) {
          1 /
            (
              -(1 + s) * m2 / (n_M - 1 - s) -
                m2 * (1 - s) / (n_M - 1 - s)
            ) *
            (s / (n_M - 1)) *
            (1 - s / (n_M - 1)) /
            (b - a + 1) *
            dgamma(
              -m2,
              shape = eta_0,
              scale = scale_val
            )
        }
        
        fm1m2_zlx <- sum(
          sapply(a:b, part10)
        )
      } else {
        part20 <- function(c) {
          1 /
            (
              -(1 + s) * c / (n_M - 1 - s) -
                m2
            ) /
            (
              -(1 + s) * c / (n_M - 1 - s) -
                c * (1 - s) / (n_M - 1 - s)
            ) *
            (
              (n_M - 1 - s) * (n_M - 2 - s) /
                ((n_M - 1) * (n_M - 2))
            ) /
            2 /
            (b - a + 1) *
            dgamma(
              -c,
              shape = eta_0,
              scale = scale_val
            )
        }
        
        part2 <- 0
        
        for (s in a:b) {
          upperbound <- min(
            0,
            -m1 * (n_M - 1 - s) / (1 + s)
          )
          lowerbound <- m2 *
            (n_M - 1 - s) /
            (1 - s)
          
          part21 <- cubature::cubintegrate(
            f = part20,
            lower = lowerbound,
            upper = upperbound,
            method = "hcubature"
          )$integral
          
          part2 <- part2 + part21
        }
        
        fm1m2_zlx <- part2
      }
    } else {
      if (m1 == m.min) {
        part12 <- function(s) {
          1 /
            (
              -(1 + s) * m1 / (n_M - 1 - s) -
                m1 * (1 - s) / (n_M - 1 - s)
            ) *
            (s / (n_M - 1)) *
            (1 - s / (n_M - 1)) /
            (b - a + 1) *
            dgamma(
              -m1,
              shape = eta_0,
              scale = scale_val
            )
        }
        
        fm1m2_zlx <- sum(
          sapply(a:b, part12)
        )
      } else {
        part22 <- function(c) {
          1 /
            (
              -(1 + s) * c / (n_M - 1 - s) -
                m1
            ) /
            (
              -(1 + s) * c / (n_M - 1 - s) -
                c * (1 - s) / (n_M - 1 - s)
            ) *
            (
              (n_M - 1 - s) * (n_M - 2 - s) /
                ((n_M - 1) * (n_M - 2))
            ) /
            2 /
            (b - a + 1) *
            dgamma(
              -c,
              shape = eta_0,
              scale = scale_val
            )
        }
        
        part2 <- 0
        
        for (s in a:b) {
          upperbound <- min(
            0,
            -m2 * (n_M - 1 - s) / (1 + s)
          )
          lowerbound <- m1 *
            (n_M - 1 - s) /
            (1 - s)
          
          part21 <- cubature::cubintegrate(
            f = part22,
            lower = lowerbound,
            upper = upperbound,
            method = "hcubature"
          )$integral
          
          part2 <- part2 + part21
        }
        
        fm1m2_zlx <- part2
      }
    }
    
    return(fm1m2_zlx)
  }
  
  # Calculate f(clrM | Z = 1, L, X)
  data1 <- as.matrix(
    cbind(
      clrM[, 1:2],
      rep(1, n_true),
      L,
      X,
      clrM.min
    )
  )
  f_clrm_z1l <- apply(data1, 1, f_clrm1m2_zl)
  
  # Calculate f(clrM | Z = 1, L = 0, X)
  data2 <- as.matrix(
    cbind(
      clrM[, 1:2],
      rep(1, n_true),
      rep(0, n_true),
      X,
      clrM.min
    )
  )
  f_clrm_z1l0 <- apply(data2, 1, f_clrm1m2_zl)
  
  # Calculate f(clrM | Z = 1, L = 1, X)
  data3 <- as.matrix(
    cbind(
      clrM[, 1:2],
      rep(1, n_true),
      rep(1, n_true),
      X,
      clrM.min
    )
  )
  f_clrm_z1l1 <- apply(data3, 1, f_clrm1m2_zl)
  
  # Calculate f(clrM | Z = 0, L, X)
  data4 <- as.matrix(
    cbind(
      clrM[, 1:2],
      rep(0, n_true),
      L,
      X,
      clrM.min
    )
  )
  f_clrm_z0l <- apply(data4, 1, f_clrm1m2_zl)
  
  # Calculate f(clrM | Z = 0, L = 0, X)
  data5 <- as.matrix(
    cbind(
      clrM[, 1:2],
      rep(0, n_true),
      rep(0, n_true),
      X,
      clrM.min
    )
  )
  f_clrm_z0l0 <- apply(data5, 1, f_clrm1m2_zl)
  
  # Calculate f(clrM | Z = 0, L = 1, X)
  data6 <- as.matrix(
    cbind(
      clrM[, 1:2],
      rep(0, n_true),
      rep(1, n_true),
      X,
      clrM.min
    )
  )
  f_clrm_z0l1 <- apply(data6, 1, f_clrm1m2_zl)
  
  # Calculate P(L | Z = 1, clrM, X)
  p_l0_z1 <- 1 - expit(
    gamma_0 +
      gamma_z +
      X %*% gamma_x
  )
  p_l1_z1 <- 1 - p_l0_z1
  
  p_l_z1clrm <- f_clrm_z1l * p_l_z1 /
    (
      f_clrm_z1l0 * p_l0_z1 +
        f_clrm_z1l1 * p_l1_z1
    )
  
  # Calculate P(L | Z = 0, clrM, X)
  p_l0_z0 <- 1 - expit(
    gamma_0 +
      X %*% gamma_x
  )
  p_l1_z0 <- 1 - p_l0_z0
  
  p_l_z0clrm <- f_clrm_z0l * p_l_z0 /
    (
      f_clrm_z0l0 * p_l0_z0 +
        f_clrm_z0l1 * p_l1_z0
    )
  
  # Estimate eta_3
  p_z <- expit(
    alpha_0 +
      X %*% alpha_x
  )
  
  w_1 <- p_l_z1 /
    (p_z * p_l_z1clrm)
  
  eta_3 <- mean(
    Z * Y * w_1
  )
  
  # Estimate eta_1
  E_Y_z1 <- expit(
    beta_0 +
      beta_z +
      beta_l * L +
      X %*% beta_x +
      clrM %*% beta_M
  )
  
  E_Y_z0 <- expit(
    beta_0 +
      beta_l * L +
      X %*% beta_x +
      clrM %*% beta_M
  )
  
  w_2 <- E_Y_z1 * p_l_z1 /
    (
      (1 - p_z) *
        E_Y_z0 *
        p_l_z0clrm
    )
  
  eta_1 <- mean(
    (1 - Z) * Y * w_2
  )
  
  # Calculate the IIE
  IIE <- eta_3 - eta_1
  
  return(IIE)
}
