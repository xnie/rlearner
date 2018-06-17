#' X-learner, as proposed by Künzel, Sekhon, Bickel, and Yu 2017, implemented via xgboost (gradient boosting)
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the treatment variable (0 or 1)
#' @param nfolds.1 number of folds for learning E[Y|X,W=1]
#' @param nfolds.0 number of folds for learning E[Y|X,W=0]
#' @param nfolds.W number of folds for learning E[W|X]
#' @param y.1.pred pre-computed estimates on E[Y|X,W=1] corresponding to the input X. xboost will compute it internally if not provided.
#' @param y.0.pred pre-computed estimates on E[Y|X,W=0] corresponding to the input X. xboost will compute it internally if not provided.
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
#' xboost.fit = xboost(X, Y, W)
#' xboost.est = predict(xboost.fit, X)
#' }
#'
#' @export
xboost = function(X, Y, W,
                  nfolds.1=NULL,
                  nfolds.0=NULL,
                  nfolds.W=NULL,
                  y.1.pred=NULL,
                  y.0.pred=NULL,
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

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(nfolds.1)) {
    nfolds.1 = floor(max(3, min(10,nobs.1/4)))
  }

  if (is.null(nfolds.0)) {
    nfolds.0 = floor(max(3, min(10,nobs.0/4)))
  }

  if (is.null(nfolds.W)) {
    nfolds.W = floor(max(3, min(10,nobs/4)))
  }

  if (is.null(y.1.pred)){
    t.1.fit = cvboost(X.1, Y.1, objective="reg:linear", nfolds = nfolds.1, nthread=nthread)
    y.1.pred = predict(t.1.fit, newx=X)
  }

  if (is.null(y.0.pred)){
    t.0.fit = cvboost(X.0, Y.0, objective="reg:linear", nfolds = nfolds.0, nthread=nthread)
    y.0.pred = predict(t.0.fit, newx=X)
  }

  D.1 = Y.1 - y.0.pred[W==1]
  D.0 = y.1.pred[W==0] - Y.0

  x.1.fit = cvboost(X.1, D.1, objective="reg:linear", nfolds = nfolds.1, nthread=nthread)
  x.0.fit = cvboost(X.0, D.0, objective="reg:linear", nfolds = nfolds.0, nthread=nthread)

  tau.1.pred = predict(x.1.fit, newx=X)
  tau.0.pred = predict(x.0.fit, newx=X)

  w.fit = cvboost(X, W, objective="binary:logistic", nfolds = nfolds.W, nthread=nthread)
  w.hat = predict(w.fit)

  tau.hat = tau.1.pred * (1-w.hat) + tau.0.pred * w.hat

  ret = list(t.1.fit = t.1.fit,
             t.0.fit = t.0.fit,
             x.1.fit = x.1.fit,
             x.0.fit = x.0.fit,
             w.fit = w.fit,
             y.1.pred = y.1.pred,
             y.0.pred = y.0.pred,
             tau.1.pred = tau.1.pred,
             tau.0.pred = tau.0.pred,
             w.hat = w.hat,
             tau.hat = tau.hat)
  class(ret) <- "xboost"
  ret

}

#' predict for xboost
#'
#' get estimated tau(x) using the trained xboost model
#'
#' @param object a xboost object
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
#' xboost.fit = xboost(X, Y, W)
#' xboost.est = predict(xboost.fit, X)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.xboost <- function(object,
                           newx=NULL,
                           ...) {
  if (!is.null(newx)) {
    tau.1.pred = predict(object$x.1.fit, newx=newx)
    tau.0.pred = predict(object$x.0.fit, newx=newx)
    w.hat = predict(object$w.fit, newx=newx)
    tau.hat = tau.1.pred * (1-w.hat) + tau.0.pred * w.hat
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
