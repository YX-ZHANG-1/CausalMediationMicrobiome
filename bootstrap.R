###############################################################################################################
# Estimate the IIE based on bootstrap samples

# Input: 
# data: a monte carlo sample
# nbts: number of bootstrap times 
# n_M: number of potential mediators
# penalty: type of penalty function
# tune: tunning parameter selection method
# iter1: maximum iteration times for iterations
# threshold: the threshold value for the estimated values


# Output:
# result_tbl: estimated values of IIE based on bootstrap samples
################################################################################################################



IPW.bootstrap <- function(data, nbts, n_M, penalty, tune, iter1, threshold){
  

  n_subj <- dim(data)[1]
  
  IIE_bs_all <- rep(NA, nbts)
 
    
    for (bs in 1:nbts){
  
      condition1 <- TRUE
      condition2 <- TRUE
      condition3 <- TRUE
      counter <- 0
    
    
      while((condition1| condition2| condition3) & (counter < iter1)){
      #resample the data with replacement
      seed <- (bs - 1) * iter1 + counter
      set.seed(seed)
      temp_df <- data %>% sample_n(size = n_subj, replace = TRUE)
      
      #count the number of iterations
      counter <- counter + 1
    
      #estimate the Y model 
      x = as.matrix(temp_df[,-1])
      y = temp_df[,1]
      
      optobj_M = rep(Inf, n_M)
      coef_M <- matrix(0, (n_M + 5), n_M)
      
      for (j in 1:n_M){
        #for each j, estimate the Y model by adding the constraint that the estimated model contains at least L,Z,X and the jth mediator
        penalty.factor0 <- c(rep(0, 4), rep(1, n_M))
        penalty.factor0[4+j] <- 0
        
        a <- list()
        a <- tryCatch({
          model <- tune.fit(x, y, family = 'binomial', penalty = penalty, tune = tune, nfolds = 5, type.measure = "deviance", penalty.factor = penalty.factor0)
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
        IIE_hat_bs <- IIE.IPW.est(temp_df = temp_df, beta_hat = beta_hat)
        
        list(val = IIE_hat_bs, error = NA)
      },
      
      error = function(e) { 
        list(val = NA, error = e)
      })
      
      #if there is any error in the estimation of IIE, then regenerate the data
      signal2 <- b$error[1]
      condition2 <- (is.na(signal2) == FALSE)
      
      if (!condition2){
        IIE_hat_bs <- b$val
        
        #if the estimated value of IIE is NA or its absolute value is larger than the threshold value, then regenerate the data
        condition3 <- (is.na(IIE_hat_bs) == TRUE | abs(IIE_hat_bs)>threshold)
      }
      
    }
    
    
    if (counter < iter1){
      IIE_hat_bs <- b$val
      IIE_bs_all[bs] <- IIE_hat_bs 
}
    
}
  
  
  result_tbl <- list(IIE_bs_all = IIE_bs_all) %>% as_tibble()

 
return(result_tbl)

}
