#' @include utils.R
#'
#' @title U-learner implemented via xgboost (boosting)
#'
#' @description U-learner as proposed by Künzel, Sekhon, Bickel, and Yu (2017), implemented via xgboost (boosting)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds number of folds used for cross fitting and cross validation
#' @param p_hat pre-computed estimates on E[W|X] corresponding to the input x. uboost will compute it internally if not provided.
#' @param m_hat pre-computed estimates on E[Y|X] corresponding to the input x. uboost will compute it internally if not provided.
#' @param cutoff the threshold to cutoff propensity estimate
#' @param ntrees_max the maximum number of trees to grow for xgboost
#' @param num_search_rounds the number of random sampling of hyperparameter combinations for cross validating on xgboost trees
#' @param print_every_n the number of iterations (in each iteration, a tree is grown) by which the code prints out information
#' @param early_stopping_rounds the number of rounds the test error stops decreasing by which the cross validation in finding the optimal number of trees stops
#' @param nthread the number of threads to use. The default is NULL, which uses all available threads
#' @param verbose boolean; whether to print statistic
#' @param bayes_opt if set to TRUE, use bayesian optimization to do hyper-parameter search in xgboost. if set to FALSE, randomly draw combinations of hyperparameters to search from (as specified by num_search_rounds). Default is FALSE.
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' uboost_fit = uboost(x, y, w)
#' uboost_est = predict(uboost_fit, x)
#' }
#'
#' @export
uboost= function(x, w, y,
                 k_folds=NULL,
                 p_hat = NULL,
                 m_hat = NULL,
                 cutoff = 0.05,
                 ntrees_max = 1000,
                 num_search_rounds = 10,
                 print_every_n = 100,
                 early_stopping_rounds = 10,
                 nthread = NULL,
                 verbose = FALSE,
                 bayes_opt = FALSE) {

  c(x, w, y) %<-% sanitize_input(x,w,y)

  nobs = nrow(x)
  pobs = ncol(x)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(w)/4)))
  }

  if (is.null(m_hat)){
    y_fit = cvboost(x,
                    y,
                    objective = "reg:linear",
                    k_folds = k_folds,
                    ntrees_max = ntrees_max,
                    num_search_rounds = num_search_rounds,
                    print_every_n = print_every_n,
                    early_stopping_rounds = early_stopping_rounds,
                    nthread = nthread,
                    verbose = verbose,
                    bayes_opt = bayes_opt)

    m_hat = predict(y_fit)
  }
  else {
    y_fit = NULL
  }

  if (is.null(p_hat)){
    w_fit = cvboost(x,
                    w,
                    objective = "binary:logistic",
                    k_folds = k_folds,
                    ntrees_max = ntrees_max,
                    num_search_rounds = num_search_rounds,
                    print_every_n = print_every_n,
                    early_stopping_rounds = early_stopping_rounds,
                    nthread = nthread,
                    verbose = verbose,
                    bayes_opt = bayes_opt)
    p_hat = predict(w_fit)
  }
  else{
    w_fit = NULL
  }

  p_hat = pmax(cutoff, pmin(1 - cutoff, p_hat))

  y_tilde = y - m_hat
  w_tilde = w - p_hat
  pseudo_outcome = y_tilde/w_tilde

  tau_fit = cvboost(x,
                    pseudo_outcome,
                    objective = "reg:linear",
                    k_folds = k_folds,
                    ntrees_max = ntrees_max,
                    num_search_rounds = num_search_rounds,
                    print_every_n = print_every_n,
                    early_stopping_rounds = early_stopping_rounds,
                    nthread = nthread,
                    verbose = verbose,
                    bayes_opt = bayes_opt)

  ret = list(tau_fit = tau_fit,
             pseudo_outcome = pseudo_outcome,
             w_fit = w_fit,
             y_fit = y_fit,
             p_hat = p_hat,
             m_hat = m_hat)
  class(ret) <- "uboost"
  ret
}

#' predict for uboost
#'
#' get estimated tau(x) using the trained uboost model
#'
#' @param object a uboost object
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
#' uboost_fit = uboost(x, w, y)
#' uboost_est = predict(uboost_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.uboost<- function(object,
                          newx = NULL,
                          ...) {
  predict(object$tau_fit, newx = newx)
}
