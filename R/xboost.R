#' @include utils.R
#'
#' @title X-learner implemented via xgboost (boosting)
#'
#' @description X-learner as proposed by Künzel, Sekhon, Bickel, and Yu (2017), implemented via xgboost (boosting)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds_mu1 number of folds for learning E[Y|X,W=1]
#' @param k_folds_mu0 number of folds for learning E[Y|X,W=0]
#' @param k_folds_p number of folds for learning E[W|X]
#' @param mu1_hat pre-computed estimates on E[Y|X,W=1] corresponding to the input x. xboost will compute it internally if not provided
#' @param mu0_hat pre-computed estimates on E[Y|X,W=0] corresponding to the input x. xboost will compute it internally if not provided
#' @param p_hat pre-computed estimates on E[W|X] corresponding to the input x. xboost will compute it internally if not provided
#' @param ntrees_max the maximum number of trees to grow for xgboost
#' @param num_search_rounds the number of random sampling of hyperparameter combinations for cross validating on xgboost trees
#' @param print_every_n the number of iterations (in each iteration, a tree is grown) by which the code prints out information
#' @param early_stopping_rounds the number of rounds the test error stops decreasing by which the cross validation in finding the optimal number of trees stops
#' @param nthread the number of threads to use. The default is NULL, which uses all available threads
#' @param verbose boolean; whether to print statistic
#' @param bayes_opt if set to TRUE, use bayesian optimization to do hyper-parameter search in xgboost. if set to FALSE, randomly draw combinations of hyperparameters to search from (as specified by num_search_rounds). default is FALSE.
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' xboost_fit = xboost(x, w, y)
#' xboost_est = predict(xboost_fit, x)
#' }
#'
#' @export
xboost = function(x, w, y,
                  k_folds_mu1=NULL,
                  k_folds_mu0=NULL,
                  k_folds_p=NULL,
                  mu1_hat=NULL,
                  mu0_hat=NULL,
                  p_hat=NULL,
                  ntrees_max=1000,
                  num_search_rounds=10,
                  print_every_n=100,
                  early_stopping_rounds=10,
                  nthread=NULL,
                  verbose=FALSE,
                  bayes_opt=FALSE) {

  c(x, w, y) %<-% sanitize_input(x,w,y)

  x_1 = x[which(w==1),]
  x_0 = x[which(w==0),]

  y_1 = y[which(w==1)]
  y_0 = y[which(w==0)]

  nobs_1 = nrow(x_1)
  nobs_0 = nrow(x_0)

  nobs = nrow(x)
  pobs = ncol(x)

  if (is.null(k_folds_mu1)) {
    k_folds_mu1 = floor(max(3, min(10,nobs_1/4)))
  }

  if (is.null(k_folds_mu0)) {
    k_folds_mu0 = floor(max(3, min(10,nobs_0/4)))
  }

  if (is.null(k_folds_p)) {
    k_folds_p = floor(max(3, min(10,nobs/4)))
  }

  if (is.null(mu1_hat)){
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
    mu1_hat = predict(t_1_fit, newx = x)
  }

  if (is.null(mu0_hat)){
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
    mu0_hat = predict(t_0_fit, newx = x)
  }

  d_1 = y_1 - mu0_hat[w==1]
  d_0 = mu1_hat[w==0] - y_0

  x_1_fit = cvboost(x_1,
                    d_1,
                    objective="reg:linear",
                    k_folds = k_folds_mu1,
                    ntrees_max = ntrees_max,
                    num_search_rounds = num_search_rounds,
                    print_every_n = print_every_n,
                    early_stopping_rounds = early_stopping_rounds,
                    nthread = nthread,
                    verbose = verbose,
                    bayes_opt = bayes_opt)

  x_0_fit = cvboost(x_0,
                    d_0,
                    objective="reg:linear",
                    k_folds = k_folds_mu0,
                    ntrees_max = ntrees_max,
                    num_search_rounds = num_search_rounds,
                    print_every_n = print_every_n,
                    early_stopping_rounds = early_stopping_rounds,
                    nthread = nthread,
                    verbose = verbose,
                    bayes_opt = bayes_opt)

  tau_1_pred = predict(x_1_fit, newx = x)
  tau_0_pred = predict(x_0_fit, newx = x)


  w_fit = cvboost(x,
                  w,
                  objective = "binary:logistic",
                  k_folds = k_folds_p,
                  ntrees_max = ntrees_max,
                  num_search_rounds = num_search_rounds,
                  print_every_n = print_every_n,
                  early_stopping_rounds = early_stopping_rounds,
                  nthread = nthread,
                  verbose = verbose,
                  bayes_opt = bayes_opt)

  p_hat = predict(w_fit)

  tau_hat = tau_1_pred * (1 - p_hat) + tau_0_pred * p_hat

  ret = list(t_1_fit = t_1_fit,
             t_0_fit = t_0_fit,
             x_1_fit = x_1_fit,
             x_0_fit = x_0_fit,
             w_fit = w_fit,
             mu1_hat = mu1_hat,
             mu0_hat = mu0_hat,
             tau_1_pred = tau_1_pred,
             tau_0_pred = tau_0_pred,
             p_hat = p_hat,
             tau_hat = tau_hat)
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
#' xboost_fit = xboost(x, w, y)
#' xboost_est = predict(xboost_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.xboost <- function(object,
                           newx = NULL,
                           ...) {
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    tau_1_pred = predict(object$x_1_fit, newx = newx)
    tau_0_pred = predict(object$x_0_fit, newx = newx)
    p_hat = predict(object$w_fit, newx = newx)
    tau_hat = tau_1_pred * (1 - p_hat) + tau_0_pred * p_hat
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
