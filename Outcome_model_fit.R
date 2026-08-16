###############################################################################
# Estimate an outcome model with penalized regression and tune the penalty
# parameter using cross-validation or an information criterion.
#
# Main function: tune.fit
#
# Arguments:
# x:                   predictor matrix
# y:                   outcome variable
# family:              outcome distribution
# penalty:             penalty function used to estimate the outcome model
# concavity.parameter: concavity parameter for SCAD or MCP
# tune:                tuning-parameter selection method
# nfolds:              number of cross-validation folds
# type.measure:         loss measure used for cross-validation
# gamma.ebic:           gamma parameter used by the extended BIC
# penalty.factor:       penalty multiplier for each coefficient
# max.iter:             maximum number of iterations
#
# Value:
# A list containing:
# coef.beta: estimated regression coefficients
# optobj:    optimal value of the tuning objective
# lambda:    selected tuning parameter
###############################################################################


# Replicate a matrix in the row and column directions
repmat <- function(X, m, n) {
  X <- as.matrix(X)
  
  mx <- nrow(X)
  nx <- ncol(X)
  
  matrix(
    t(matrix(X, mx, nx * n)),
    nrow = mx * m,
    ncol = nx * n,
    byrow = TRUE
  )
}


# Calculate the loss for each candidate coefficient vector
loglik <- function(X, y, beta, family) {
  K <- ncol(beta)
  link <- cbind(1, X) %*% beta
  yrep <- repmat(y, 1, K)
  
  if (family == "gaussian") {
    return(apply((yrep - link)^2, 2, sum))
  }
  
  if (family == "poisson") {
    return(apply(exp(link) - yrep * link, 2, sum))
  }
  
  if (family == "binomial") {
    return(apply(log(1 + exp(link)) - yrep * link, 2, sum))
  }
}


# Calculate the number of nonzero coefficients
getdf <- function(coef.beta) {
  apply(abs(coef.beta) > 1e-10, 2, sum)
}


# Fit and tune the penalized outcome model
tune.fit <- function(
    x,
    y,
    family = c("gaussian", "binomial", "poisson", "cox"),
    penalty = c("SCAD", "MCP", "lasso"),
    concavity.parameter = switch(penalty, SCAD = 3.7, 3),
    tune = c("cv", "aic", "bic", "ebic"),
    nfolds = 10,
    type.measure = c("deviance", "class", "auc", "mse", "mae"),
    gamma.ebic = 1,
    penalty.factor,
    max.iter = 10000
) {
  if (is.null(x) || is.null(y)) {
    stop("The data are missing.")
  }
  
  this.call <- match.call()
  
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  tune <- match.arg(tune)
  type.measure <- match.arg(type.measure)
  
  if (!is.numeric(concavity.parameter)) {
    stop("concavity.parameter must be numeric.")
  }
  
  if (!is.numeric(nfolds)) {
    stop("nfolds must be numeric.")
  }
  
  if (tune == "cv") {
    # Select lambda using cross-validation
    if (penalty == "lasso") {
      cv.fit <- cv.glmnet(
        x,
        y,
        family = family,
        type.measure = type.measure,
        nfolds = nfolds,
        penalty.factor = penalty.factor
      )
      
      coef.beta <- coef(cv.fit, s = "lambda.min")
      reg.fit <- cv.fit$glmnet.fit
      lambda <- cv.fit$lambda.min
      lambda.ind <- which(cv.fit$lambda == cv.fit$lambda.min)
      optobj <- min(cv.fit$cvm)
    } else if (family != "cox") {
      cv.fit <- cv.ncvreg(
        x,
        y,
        family = family,
        penalty = penalty,
        gamma = concavity.parameter,
        nfolds = nfolds,
        penalty.factor = penalty.factor,
        max.iter = max.iter
      )
      
      coef.beta <- coef(cv.fit, s = "lambda.min")
      reg.fit <- cv.fit$fit
      lambda <- cv.fit$lambda.min
      lambda.ind <- which(cv.fit$lambda == cv.fit$lambda.min)
      optobj <- min(cv.fit$cve)
    } else {
      cv.fit <- cv.ncvsurv(
        x,
        y,
        family = family,
        penalty = penalty,
        gamma = concavity.parameter,
        nfolds = nfolds,
        penalty.factor = penalty.factor,
        max.iter = max.iter
      )
      
      coef.beta <- coef(cv.fit, s = "lambda.min")
      reg.fit <- cv.fit$fit
      lambda <- cv.fit$lambda.min
      lambda.ind <- which(cv.fit$lambda == cv.fit$lambda.min)
      optobj <- min(cv.fit$cve)
    }
  } else {
    # Select lambda using an information criterion
    n <- nrow(x)
    
    if (penalty == "lasso") {
      reg.fit <- glmnet(
        x,
        y,
        family = family,
        penalty.factor = penalty.factor
      )
      
      # Extract coefficients, including the intercept, for every lambda
      coef.beta <- rbind(
        reg.fit$a0,
        as.matrix(reg.fit$beta)
      )
      
      dev <- deviance(reg.fit)
      reg.df <- reg.fit$df
    } else if (family != "cox") {
      reg.fit <- ncvreg(
        x,
        y,
        family = family,
        penalty = penalty,
        gamma = concavity.parameter,
        penalty.factor = penalty.factor,
        max.iter = max.iter,
        trace = FALSE
      )
      
      # Extract coefficients, including the intercept, for every lambda
      coef.beta <- reg.fit$beta
      dev <- loglik(x, y, coef.beta, family = family)
      reg.df <- getdf(coef.beta[-1, , drop = FALSE])
    } else {
      reg.fit <- ncvsurv(
        x,
        y,
        family = family,
        penalty = penalty,
        gamma = concavity.parameter,
        penalty.factor = penalty.factor,
        max.iter = max.iter
      )
      
      coef.beta <- reg.fit$beta
      dev <- 2 * reg.fit$loss
      reg.df <- getdf(coef.beta)
    }
    
    if (tune == "aic") {
      obj <- dev + 2 * reg.df
    }
    
    if (tune == "bic") {
      obj <- dev + log(n) * reg.df
    }
    
    if (tune == "ebic") {
      obj <- dev +
        log(n) * reg.df +
        2 * gamma.ebic * log(choose(ncol(x), reg.df))
    }
    
    lambda.ind <- which.min(obj)
    coef.beta <- coef.beta[, lambda.ind]
    lambda <- reg.fit$lambda[lambda.ind]
    optobj <- min(obj)
  }
  
  return(
    list(
      coef.beta = coef.beta,
      optobj = optobj,
      lambda = lambda
    )
  )
}
