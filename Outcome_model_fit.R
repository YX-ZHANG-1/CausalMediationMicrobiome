###############################################################################################################
# Estimate outcome model

# Main function: tune.fit
# Input: 
# x: covariates 
# y: outcome variable 
# family: distribution family
# penalty: type of penalty function
# concavity.parameter: concavity parameter 
# tune: tunning parameter selection method
# nfolds: the number of folds for cross-validation
# type.measure: loss to use for cross-validation
# gamma.ebic: choice of gamma for extended BIC method
# penalty.factor: a number that multiplies lambda to allow differential shrinkage


# Output:
# coef.beta: estimated values of coefficients 
# optobj: optimal objective function
################################################################################################################

#define some subfunctions which will be used by other functions
repmat <- function(X, m, n) {
  X <- as.matrix(X)
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}


loglik <- function(X, y, beta, family) {
  K <- dim(beta)[2]
  link <- cbind(1, X) %*% beta
  yrep <- repmat(y, 1, K)
  if (family == "gaussian") 
    return(apply((yrep - link)^2, 2, sum))
  if (family == "poisson") 
    return(apply(exp(link) - yrep * link, 2, sum))
  if (family == "binomial") 
    return(apply(log(1 + exp(link)) - yrep * link, 2, sum))
}


getdf <- function(coef.beta) {
  apply(abs(coef.beta) > 1e-10, 2, sum)
}


#the main function for estimating the Y model
tune.fit <- function(x, y, family = c("gaussian", "binomial", "poisson", "cox"), penalty = c("SCAD", "MCP", "lasso"), concavity.parameter = switch(penalty, SCAD = 3.7, 3), tune = c("cv", "aic", "bic", "ebic"), nfolds = 10, 
                     type.measure = c("deviance", "class", "auc", "mse", "mae"), gamma.ebic = 1, penalty.factor) {
  
  if (is.null(x) || is.null(y)) 
    stop("The data is missing!")
  
  this.call = match.call()
  family = match.arg(family)
  penalty = match.arg(penalty)
  if (class(concavity.parameter) != "numeric") 
    stop("concavity.parameter must be numeric!")
  tune = match.arg(tune)
  if (class(nfolds) != "numeric") 
    stop("nfolds must be numeric!")
  type.measure = match.arg(type.measure)
  
  
  if (tune == "cv") {
    if (penalty == "lasso" ) {
      cv.fit = cv.glmnet(x, y, family = family, type.measure = type.measure, nfolds = nfolds, penalty.factor = penalty.factor)
      coef.beta = coef(cv.fit, s = "lambda.min") 
      reg.fit = cv.fit$glmnet.fit
      lambda = cv.fit$lambda.min
      lambda.ind = which(cv.fit$lambda == cv.fit$lambda.min)
      optobj = min(cv.fit$cvm)
      
    } else if (family != 'cox') {
      cv.fit = cv.ncvreg(x, y, family = family, penalty = penalty, gamma = concavity.parameter, nfolds = nfolds, penalty.factor = penalty.factor, max.iter=10000)
      coef.beta = coef(cv.fit, s = "lambda.min")
      reg.fit = cv.fit$fit
      lambda = cv.fit$lambda.min
      lambda.ind = which(cv.fit$lambda == cv.fit$lambda.min)
      optobj = min(cv.fit$cve)
    } else {
      cv.fit = cv.ncvsurv(x, y, family = family, penalty = penalty, gamma = concavity.parameter, nfolds = nfolds, penalty.factor = penalty.factor, max.iter=10000)
      coef.beta = coef(cv.fit, s = "lambda.min")
      reg.fit = cv.fit$fit
      lambda = cv.fit$lambda.min
      lambda.ind = which(cv.fit$lambda == cv.fit$lambda.min)
      optobj = min(cv.fit$cve)
    }
  } else {
    n = nrow(x)
    if (penalty == "lasso" ) {
      reg.fit = glmnet(x, y, family = family, penalty.factor = penalty.factor)
      coef.beta = rbind(reg.fit$a0,as.matrix(reg.fit$beta))  # extract coefficients at all values of lambda,  including the intercept
      dev = deviance(reg.fit)
      reg.df = reg.fit$df
    } else {
      if(family != 'cox'){
        
        reg.fit = ncvreg(x, y, family = family, penalty = penalty, gamma = concavity.parameter, penalty.factor = penalty.factor, max.iter=10000, trace=FALSE)
        coef.beta = reg.fit$beta  # extract coefficients at all values of lambda, including the intercept
        dev = loglik(x, y, coef.beta, family = family)
        reg.df = getdf(coef.beta[-1, , drop = FALSE])
      } else {
        reg.fit = ncvsurv(x, y, family = family, penalty = penalty, gamma = concavity.parameter, penalty.factor = penalty.factor, max.iter=10000)
        coef.beta = reg.fit$beta  # extract coefficients at all values of lambda, including the intercept
        dev = 2*reg.fit$loss
        reg.df = getdf(coef.beta)
      }
    }
    
    if (tune == "aic") {
      obj = dev + 2 * reg.df
    }
    if (tune == "bic") {
      obj = dev + log(n) * reg.df
    }
    if (tune == "ebic") {
      obj = dev + log(n) * reg.df + 2 * gamma.ebic * log(choose(dim(x)[2], reg.df))
    }
    lambda.ind = which.min(obj)
    coef.beta = coef.beta[, lambda.ind]
    lambda = reg.fit$lambda[lambda.ind]
    optobj = min(obj)
  }
  

  return(list(coef.beta = coef.beta, optobj = optobj, lambda = lambda))
}
