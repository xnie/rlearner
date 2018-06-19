#' S-learner, as proposed by Imai and Ratkovic 2013, implemented via xgboost (gradient boosting)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds number of folds for cross validation
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
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' sboost.fit = sboost(x, w, y)
#' sboost.est = predict(sboost.fit, x)
#' }
#'
#' @export
sboost = function(X, W, Y,
                  k_folds = NULL,
                  ntrees.max=1000,
                  num.search.rounds=10,
                  print.every.n=100,
                  early.stopping.rounds=10,
                  nthread=NULL,
                  bayes.opt=FALSE){

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,nobs/4)))
  }

  s.fit = cvboost(cbind(X, (W-0.5)*X, (W-0.5)),
                  Y,
                  objective="reg:linear",
                  k_folds=k_folds,
                  ntrees.max=ntrees.max,
                  num.search.rounds=num.search.rounds,
                  print.every.n=print.every.n,
                  early.stopping.rounds=early.stopping.rounds,
                  nthread=nthread,
                  bayes.opt=bayes.opt)


  mu0.hat = predict(s.fit, newx=cbind(X, (0-0.5)*X, (0-0.5)))
  mu1.hat = predict(s.fit, newx=cbind(X, (1-0.5)*X, (1-0.5)))
  tau.hat = mu1.hat - mu0.hat

  ret = list(s.fit = s.fit,
             mu0.hat = mu0.hat,
             mu1.hat = mu1.hat,
             tau.hat = tau.hat)

  class(ret) <- "sboost"
  ret
}

#' predict for sboost
#'
#' get estimated tau(x) using the trained sboost model
#'
#' @param object a sboost object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param ... additional arguments (currently not used)
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' sboost.fit = sboost(x, w, y)
#' sboost.est = predict(sboost.fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.sboost <- function(object,
                           newx=NULL,
                           ...) {
  if (!is.null(newx)) {
    mu0.hat = predict(object$s.fit, newx=cbind(newx, (0-0.5)*newx, (0-0.5)))
    mu1.hat = predict(object$s.fit, newx=cbind(newx, (1-0.5)*newx, (1-0.5)))
    tau.hat = mu1.hat - mu0.hat
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
