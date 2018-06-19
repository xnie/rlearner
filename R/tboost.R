#' T-learner, implemented via xgboost (gradient boosting)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds_mu1 number of folds for cross validation for the treated
#' @param k_folds_mu0 number of folds for cross validation for the control
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
#' tboost_fit = tboost(x, w, y)
#' tboost_est = predict(tboost_fit, x)
#' }
#'
#' @export

tboost = function(x, w, y,
                  alpha = 1,
                  k_folds_mu1 = NULL,
                  k_folds_mu0 = NULL,
                  ntrees_max = 1000,
                  num_search_rounds = 10,
                  print_every_n = 100,
                  early_stopping_rounds = 10,
                  nthread = NULL,
                  verbose = FALSE,
                  bayes_opt = FALSE) {

  x_1 = x[which(w==1),]
  x_0 = x[which(w==0),]

  y_1 = y[which(w==1)]
  y_0 = y[which(w==0)]

  nobs_1 = nrow(x_1)
  nobs_0 = nrow(x_0)

  pobs = ncol(x)

  if (is.null(k_folds_mu1)) {
    k_folds_mu1 = floor(max(3, min(10,nobs_1/4)))
  }

  if (is.null(k_folds_mu0)) {
    k_folds_mu0 = floor(max(3, min(10,nobs_0/4)))
  }

  t_1_fit = cvboost(x_1,
                    y_1,
                    objective="reg:linear",
                    k_folds = k_folds_mu1,
                    ntrees_max = ntrees_max,
                    num_search_rounds = num_search_rounds,
                    print_every_n = print_every_n,
                    early_stopping_rounds = early_stopping_rounds,
                    nthread = nthread,
                    verbose = verbose,
                    bayes_opt = bayes_opt)

  t_0_fit = cvboost(x_0,
                    y_0,
                    objective = "reg:linear",
                    k_folds = k_folds_mu0,
                    ntrees_max = ntrees_max,
                    num_search_rounds = num_search_rounds,
                    print_every_n = print_every_n,
                    early_stopping_rounds = early_stopping_rounds,
                    nthread = nthread,
                    verbose = verbose,
                    bayes_opt = bayes_opt)

  y_1_pred = predict(t_1_fit, newx=x)
  y_0_pred = predict(t_0_fit, newx=x)

  tau_hat = y_1_pred - y_0_pred

  ret = list(t_1_fit = t_1_fit,
             t_0_fit = t_0_fit,
             y_1_pred = y_1_pred,
             y_0_pred = y_0_pred,
             tau_hat = tau_hat)
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
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' tboost_fit = tboost(x, w, y)
#' tboost_est = predict(tboost_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.tboost <- function(object,
                           newx = NULL,
                           ...) {
  if (!is.null(newx)) {
    y_1_pred = predict(object$t_1_fit, newx = newx)
    y_0_pred = predict(object$t_0_fit, newx = newx)
    tau_hat = y_1_pred - y_0_pred
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
