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
#' @export rlassocf
rlassocf= function(X, Y, W,
                  alpha = 1,
                  nfolds=NULL,
                  lambda.choice=c("lambda.1se", "lambda.min"),
                  standardize = TRUE) {
  
  lambda.choice = match.arg(lambda.choice)
  
  nobs = nrow(X)
  pobs = ncol(X)
  
  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,length(W)/4)))
  }
  
  # fold ID for cross-validation; balance treatment assignments
  foldid = sample(rep(seq(nfolds), length = length(W)))
  
  y.fit = crossfit.cv.glmnet(X, Y, foldid = foldid, lambda.choice = lambda.choice, alpha = alpha, standardize = standardize)
  y.hat = crossfit.predict(y.fit)
  
  w.fit = crossfit.cv.glmnet(X, W, foldid = foldid, lambda.choice = lambda.choice, alpha = alpha, standardize = standardize,family="binomial", type.measure = "auc")
  w.hat = crossfit.predict(w.fit)
  
  Y.tilde = Y - y.hat
  X.tilde = cbind(as.numeric(W - w.hat) * cbind(1, X))
  
  tau.fit = crossfit.cv.glmnet(X.tilde, Y.tilde, foldid = foldid,
                               lambda.choice = lambda.choice, alpha = alpha,
                               penalty.factor = c(0, rep(1, pobs)),
                               standardize = standardize)
  tau.hat = crossfit.predict(tau.fit, cbind(1, X))
  tau.beta = coef(tau.fit)[1 + 1:(1+pobs)]
  
  return(list(tau.hat = tau.hat, y.hat = y.hat, w.hat = w.hat, tau.beta = tau.beta))
}

crossfit.predict = function(lasso.obj, X = lasso.obj$x) {
  sapply(1:length(lasso.obj$foldid), function(i) {
    sum(lasso.obj$cv.betas[lasso.obj$foldid[i],] * c(1, X[i,]))
  })
}
