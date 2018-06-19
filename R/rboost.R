#' R-learner, as proposed by Nie and Wager 2017, implemented via xgboost (gradient boosting)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds number of folds used for cross fitting and cross validation
#' @param p_hat pre-computed estimates on E[W|X] corresponding to the input X. rboost will compute it internally if not provided.
#' @param m_hat pre-computed estimates on E[Y|X] corresponding to the input X. rboost will compute it internally if not provided.
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
#' rboost.fit = rboost(x, w, y)
#' rboost.est = predict(rboost.fit, x)
#' }
#'
#' @export
rboost= function(X, W, Y,
                 k_folds=NULL,
                 p_hat = NULL,
                 m_hat = NULL,
                 ntrees.max=1000,
                 num.search.rounds=10,
                 print.every.n=100,
                 early.stopping.rounds=10,
                 nthread=NULL,
                 bayes.opt=FALSE) {

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(W)/4)))
  }

  if (is.null(m_hat)){
    y.fit = cvboost(X,
                    Y,
                    objective="reg:linear",
                    k_folds=k_folds,
                    ntrees.max=ntrees.max,
                    num.search.rounds=num.search.rounds,
                    print.every.n=print.every.n,
                    early.stopping.rounds=early.stopping.rounds,
                    nthread=nthread,
                    bayes.opt=bayes.opt)

    m_hat = predict(y.fit)
  }
  else {
    y.fit = NULL
  }

  if (is.null(p_hat)){
    w.fit = cvboost(X,
                    W,
                    objective="binary:logistic",
                    k_folds=k_folds,
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

  Y.tilde = Y - m_hat
  W.tilde = W - p_hat
  pseudo.outcome = Y.tilde/W.tilde

  weights = W.tilde^2

  tau.fit = cvboost(X,
                    pseudo.outcome,
                    objective="reg:linear",
                    weights=weights,
                    k_folds=k_folds,
                    ntrees.max=ntrees.max,
                    num.search.rounds=num.search.rounds,
                    print.every.n=print.every.n,
                    early.stopping.rounds=early.stopping.rounds,
                    nthread=nthread,
                    bayes.opt=bayes.opt)

  ret = list(tau.fit = tau.fit,
             pseudo.outcome = pseudo.outcome,
             weights = weights,
             w.fit = w.fit,
             y.fit = y.fit,
             p_hat = p_hat,
             m_hat = m_hat)
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
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rboost.fit = rboost(x, w, y)
#' rboost.est = predict(rboost.fit, x)
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
