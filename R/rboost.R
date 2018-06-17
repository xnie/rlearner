#' R-learner, as proposed by Nie and Wager 2017, implemented via xgboost (gradient boosting)
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the treatment variable (0 or 1)
#' @param nfolds number of folds used for cross fitting and cross validation
#' @param w.hat pre-computed estimates on E[W|X] corresponding to the input X. rboost will compute it internally if not provided.
#' @param y.hat pre-computed estimates on E[Y|X] corresponding to the input X. rboost will compute it internally if not provided.
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
#' rboost.fit = rboost(X, Y, W)
#' rboost.est = predict(rboost.fit, X)
#' }
#'
#' @export
rboost= function(X, Y, W,
                 nfolds=NULL,
                 w.hat = NULL,
                 y.hat = NULL,
                 ntrees.max=1000,
                 num.search.rounds=10,
                 print.every.n=100,
                 early.stopping.rounds=10,
                 nthread=NULL,
                 bayes.opt=FALSE) {

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,length(W)/4)))
  }

  if (is.null(y.hat)){
    y.fit = cvboost(X, Y, objective="reg:linear", nfolds=nfolds, nthread=nthread)
    y.hat = predict(y.fit)
  }
  else {
    y.fit = NULL
  }

  if (is.null(w.hat)){
    w.fit = cvboost(X, W, objective="binary:logistic", nfolds=nfolds, nthread=nthread)
    w.hat = predict(w.fit)
  }
  else{
    w.fit = NULL
  }

  Y.tilde = Y - y.hat
  W.tilde = W - w.hat
  pseudo.outcome = Y.tilde/W.tilde
  tau.const = NULL

  weights = W.tilde^2

  tau.fit = cvboost(X, pseudo.outcome, objective="reg:linear", weights=weights, nfolds=nfolds, nthread=nthread)

  ret = list(tau.fit = tau.fit,
             w.fit = w.fit,
             y.fit = y.fit,
             w.hat = w.hat,
             y.hat = y.hat,
             tau.const = tau.const)
  class(ret) <- "rboost"
  ret
}

#' predict for rboost
#'
#' get estimated tau(x) using the trained rboost model
#'
#' @param object a rboost object
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
#' rboost.fit = rboost(X, Y, W)
#' rboost.est = predict(rboost.fit, X)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.rboost<- function(object,
                           newx=NULL,
                           ...) {
  predict(object$tau.fit, newx=newx)
}
