#' T-learner, implemented via xgboost (gradient boosting)
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the treatment variable (0 or 1)
#' @param nfolds.1 number of folds for cross validation for the treated
#' @param nfolds.0 number of folds for cross validation for the control
#' @param ntrees.max the maximum number of trees to grow for xgboost
#' @param num.search.rounds the number of random sampling of hyperparameter combinations for cross validating on xgboost trees
#' @param print.every.n the number of iterations (in each iteration, a tree is grown) by which the code prints out information
#' @param early.stopping.rounds the number of rounds the test error stops decreasing by which the cross validation in finding the optimal number of trees stops
#' @param nthread the number of threads to use. The default is NULL, which uses all available threads
#' @param bayes.opt if set to TRUE, use bayesian optimization to do hyper-parameter search in xgboost. if set to FALSE, randomly draw combinations of hyperparameters to search from (as specified by num.search.rounds). Default is FALSE.
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#'
#' tboost.fit = tboost(X, Y, W)
#' tboost.est = predict(tboost.fit, X)
#' }
#'
#' @export

tboost = function(X, Y, W,
                  alpha = 1,
                  nfolds.1=NULL,
                  nfolds.0=NULL,
                  ntrees.max=1000,
                  num.search.rounds=10,
                  print.every.n=100,
                  early.stopping.rounds=10,
                  nthread=NULL,
                  bayes.opt=FALSE) {

  X.1 = X[which(W==1),]
  X.0 = X[which(W==0),]

  Y.1 = Y[which(W==1)]
  Y.0 = Y[which(W==0)]

  nobs.1 = nrow(X.1)
  nobs.0 = nrow(X.0)

  pobs = ncol(X)

  if (is.null(nfolds.1)) {
    nfolds.1 = floor(max(3, min(10,nobs.1/4)))
  }

  if (is.null(nfolds.0)) {
    nfolds.0 = floor(max(3, min(10,nobs.0/4)))
  }

  t.1.fit = cvboost(X.1, Y.1, objective="reg:linear", nfolds = nfolds.1, nthread=nthread)
  t.0.fit = cvboost(X.0, Y.0, objective="reg:linear", nfolds = nfolds.0, nthread=nthread)

  y.1.pred = predict(t.1.fit, newx=X)
  y.0.pred = predict(t.0.fit, newx=X)

  tau.hat = y.1.pred - y.0.pred

  ret = list(t.1.fit = t.1.fit,
             t.0.fit = t.0.fit,
             y.1.pred = y.1.pred,
             y.0.pred = y.0.pred,
             tau.hat = tau.hat)
  class(ret) <- "tboost"
  ret
}

#' predict for tboost
#'
#' get estimated tau(x) using the trained tboost model
#'
#' @param object a tboost object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param ... additional arguments (currently not used)
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#'
#' tboost.fit = tboost(X, Y, W)
#' tboost.est = predict(tboost.fit, X)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.tboost <- function(object,
                           newx=NULL,
                           ...) {
  if (!is.null(newx)) {
    y.1.pred = predict(object$t.1.fit, newx=newx)
    y.0.pred = predict(object$t.0.fit, newx=newx)
    tau.hat = y.1.pred - y.0.pred
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
