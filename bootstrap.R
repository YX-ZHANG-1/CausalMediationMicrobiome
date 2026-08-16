###############################################################################
# Estimate the interventional indirect effect (IIE) using bootstrap samples.
#
# Arguments:
# data:      Monte Carlo sample containing the analysis variables
# nbts:      number of bootstrap samples
# n_M:       number of mediator components
# penalty:   penalty function used to estimate the outcome model
# tune:      tuning-parameter selection method
# iter1:     maximum number of estimation attempts per bootstrap sample
# threshold: maximum acceptable absolute value of the IIE estimate
#
# Value:
# A tibble containing the estimated IIE for each bootstrap sample.
###############################################################################


IPW.bootstrap <- function(
    data,
    nbts,
    n_M,
    penalty,
    tune,
    iter1,
    threshold
) {
  n_subj <- nrow(data)
  IIE_bs_all <- rep(NA_real_, nbts)
  
  for (bs in seq_len(nbts)) {
    condition <- TRUE
    counter <- 0
    success_b <- FALSE
    
    IIE_attempt <- list(
      success = FALSE,
      val = NULL
    )
    
    
    while ((!success_b || condition) && counter < iter1) {
      beta_hat <- NULL
      
      # Resample the data with replacement
      seed <- (bs - 1) * iter1 + counter
      set.seed(seed)
      
      temp_df <- data %>%
        sample_n(
          size = n_subj,
          replace = TRUE
        )
      
      # Count the number of estimation attempts
      counter <- counter + 1
      
      # Prepare the outcome and predictors
      x <- as.matrix(temp_df[, -1])
      y <- temp_df[, 1]
      
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
          
          beta_hat <- as.matrix(model$coef.beta)
          
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
          IIE_hat_bs <- IIE.IPW.est(
            temp_df = temp_df,
            beta_hat = beta_hat
          )
          
          list(
            success = TRUE,
            val = IIE_hat_bs
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
        IIE_hat_bs <- IIE_attempt$val
        
        # Regenerate the data if the estimate is missing or exceeds the
        # prespecified threshold.
        condition <- is.na(IIE_hat_bs) ||
          abs(IIE_hat_bs) > threshold
      } else {
        success_b <- FALSE
        condition <- TRUE
      }
    }
    
    if (success_b && !condition) {
      IIE_bs_all[bs] <- IIE_attempt$val
    }
  }
  
  result_tbl <- tibble(
    IIE_bs_all = IIE_bs_all
  )
  
  return(result_tbl)
}
