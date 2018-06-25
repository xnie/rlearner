#' @include utils.R
#'
#' @title S-learner, implemented via xgboost (boosting)
#'
#' @description  S-learner, as proposed by Imai and Ratkovic (2013), implemented via xgboost (boosting)
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds number of folds for cross validation
#' @param ntrees_max the maximum number of trees to grow for xgboost
#' @param num_search_rounds the number of random sampling of hyperparameter combinations for cross validating on xgboost trees
#' @param print_every_n the number of iterations (in each iteration, a tree is grown) by which the code prints out information
#' @param early_stopping_rounds the number of rounds the test error stops decreasing by which the cross validation in finding the optimal number of trees stops
#' @param nthread the number of threads to use. The default is NULL, which uses all available threads
#' @param verbose boolean; whether to print statistic
#' @param bayes_opt if set to TRUE, use bayesian optimization to do hyper-parameter search in xgboost. if set to FALSE, randomly draw combinations of hyperparameters to search from (as specified by num_search_rounds). Default is FALSE. CAUTION: current implementation is in beta and is not recommended for usage yet.
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' sboost_fit = sboost(x, w, y)
#' sboost_est = predict(sboost_fit, x)
#' }
#'
#' @export
sboost = function(x, w, y,
                  k_folds = NULL,
                  ntrees_max = 1000,
                  num_search_rounds = 10,
                  print_every_n = 100,
                  early_stopping_rounds = 10,
                  nthread = NULL,
                  verbose = FALSE,
                  bayes_opt = FALSE){
  c(x, w, y) %<-% sanitize_input(x,w,y)

  nobs = nrow(x)
  pobs = ncol(x)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10, nobs/4)))
  }

  s_fit = cvboost(cbind(x, (w-0.5)*x, (w-0.5)),
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


  mu0_hat = predict(s_fit, newx = cbind(x, (0 - 0.5) * x, (0 - 0.5)))
  mu1_hat = predict(s_fit, newx = cbind(x, (1 - 0.5) * x, (1 - 0.5)))
  tau_hat = mu1_hat - mu0_hat

  ret = list(s_fit = s_fit,
             mu0_hat = mu0_hat,
             mu1_hat = mu1_hat,
             tau_hat = tau_hat)

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
#' sboost_fit = sboost(x, w, y)
#' sboost_est = predict(sboost_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.sboost <- function(object,
                           newx = NULL,
                           ...) {
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    mu0_hat = predict(object$s_fit, newx = cbind(newx, (0 - 0.5) * newx, (0 - 0.5)))
    mu1_hat = predict(object$s_fit, newx = cbind(newx, (1 - 0.5) * newx, (1 - 0.5)))
    tau_hat = mu1_hat - mu0_hat
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
