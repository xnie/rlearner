#' Estimate treatment effect via the R-lasso, as proposed by Nie and Wager (2017)
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the effect variable (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param nfolds number of folds for cross-fitting
#' @param lambda.choice how to cross-validated
#' @param standardize whether X should be rescaled before running the lasso
#'
#' @export slasso
slasso = function(X, Y, W,
                  alpha = 1,
                  nfolds=NULL,
                  lambda.choice=c("lambda.1se", "lambda.min"),
                  standardize = TRUE) {
  
  lambda.choice = match.arg(lambda.choice)
  
  nobs = nrow(X)
  pobs = ncol(X)
  
  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,nobs/4)))
  }
  
  # fold ID for cross-validation; balance treatment assignments
  foldid = sample(rep(seq(nfolds), length = nobs))
  
  X.tilde = cbind(as.numeric(2 * W - 1) * cbind(1, X), X)
  
  s.fit = cv.glmnet(X.tilde, Y, foldid = foldid, 
                      penalty.factor = c(0, rep(1, 2 * pobs)), 
                      standardize=standardize)
  
  tau.hat = 2 * cbind(1, X, 0 * X) %*% as.vector(t(coef(s.fit, s=lambda.choice)[-1]))
  
  return(list(tau.hat = tau.hat))
}