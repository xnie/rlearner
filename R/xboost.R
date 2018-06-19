#' X-learner, as proposed by Künzel, Sekhon, Bickel, and Yu 2017, implemented via xgboost (gradient boosting)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds_mu1 number of folds for learning E[Y|X,W=1]
#' @param k_folds_mu0 number of folds for learning E[Y|X,W=0]
#' @param k_folds_p number of folds for learning E[W|X]
#' @param mu1_hat pre-computed estimates on E[Y|X,W=1] corresponding to the input X. xboost will compute it internally if not provided
#' @param mu0_hat pre-computed estimates on E[Y|X,W=0] corresponding to the input X. xboost will compute it internally if not provided
#' @param p_hat pre-computed estimates on E[W|X] corresponding to the input X. xboost will compute it internally if not provided
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
#' xboost.fit = xboost(x, w, y)
#' xboost.est = predict(xboost.fit, x)
#' }
#'
#' @export
xboost = function(X, W, Y,
                  k_folds_mu1=NULL,
                  k_folds_mu0=NULL,
                  k_folds_p=NULL,
                  mu1_hat=NULL,
                  mu0_hat=NULL,
                  p_hat=NULL,
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

  if (is.null(k_folds_mu1)) {
    k_folds_mu1 = floor(max(3, min(10,nobs.1/4)))
  }

  if (is.null(k_folds_mu0)) {
    k_folds_mu0 = floor(max(3, min(10,nobs.0/4)))
  }

  if (is.null(k_folds_p)) {
    k_folds_p = floor(max(3, min(10,nobs/4)))
  }

  if (is.null(mu1_hat)){
    t.1.fit = cvboost(X.1,
                      Y.1,
                      objective="reg:linear",
                      nfolds = k_folds_mu1,
                      ntrees.max=ntrees.max,
                      num.search.rounds=num.search.rounds,
                      print.every.n=print.every.n,
                      early.stopping.rounds=early.stopping.rounds,
                      nthread=nthread,
                      bayes.opt=bayes.opt)
    mu1_hat = predict(t.1.fit, newx=X)
  }

  if (is.null(mu0_hat)){
    t.0.fit = cvboost(X.0,
                      Y.0,
                      objective="reg:linear",
                      nfolds = k_folds_mu0,
                      ntrees.max=ntrees.max,
                      num.search.rounds=num.search.rounds,
                      print.every.n=print.every.n,
                      early.stopping.rounds=early.stopping.rounds,
                      nthread=nthread,
                      bayes.opt=bayes.opt)
    mu0_hat = predict(t.0.fit, newx=X)
  }

  D.1 = Y.1 - mu0_hat[W==1]
  D.0 = mu1_hat[W==0] - Y.0

  x.1.fit = cvboost(X.1,
                    D.1,
                    objective="reg:linear",
                    nfolds = k_folds_mu1,
                    ntrees.max=ntrees.max,
                    num.search.rounds=num.search.rounds,
                    print.every.n=print.every.n,
                    early.stopping.rounds=early.stopping.rounds,
                    nthread=nthread,
                    bayes.opt=bayes.opt)

  x.0.fit = cvboost(X.0,
                    D.0,
                    objective="reg:linear",
                    nfolds = k_folds_mu0,
                    ntrees.max=ntrees.max,
                    num.search.rounds=num.search.rounds,
                    print.every.n=print.every.n,
                    early.stopping.rounds=early.stopping.rounds,
                    nthread=nthread,
                    bayes.opt=bayes.opt)

  tau.1.pred = predict(x.1.fit, newx=X)
  tau.0.pred = predict(x.0.fit, newx=X)


  if (is.null(p_hat)){
    w.fit = cvboost(X,
                    W,
                    objective="binary:logistic",
                    nfolds=k_folds_p,
                    ntrees.max=ntrees.max,
                    num.search.rounds=num.search.rounds,
                    print.every.n=print.every.n,
                    early.stopping.rounds=early.stopping.rounds,
                    nthread=nthread,
                    bayes.opt=bayes.opt)
    p_hat = predict(w.fit)
  }
  else{
    w.fit = NULL
  }

  tau.hat = tau.1.pred * (1-p_hat) + tau.0.pred * p_hat

  ret = list(t.1.fit = t.1.fit,
             t.0.fit = t.0.fit,
             x.1.fit = x.1.fit,
             x.0.fit = x.0.fit,
             w.fit = w.fit,
             mu1_hat = mu1_hat,
             mu0_hat = mu0_hat,
             tau.1.pred = tau.1.pred,
             tau.0.pred = tau.0.pred,
             p_hat = p_hat,
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
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' xboost.fit = xboost(x, w, y)
#' xboost.est = predict(xboost.fit, x)
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
    p_hat = predict(object$w.fit, newx=newx)
    tau.hat = tau.1.pred * (1-p_hat) + tau.0.pred * p_hat
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
