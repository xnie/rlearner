#' @include utils.R
#'
#' @title R-learner, implemented via xgboost (boosting)
#'
#' @description  R-learner, as proposed by Nie and Wager (2017), implemented via xgboost (boosting)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds number of folds used for cross fitting and cross validation
#' @param p_hat pre-computed estimates on E[W|X] corresponding to the input x. rboost will compute it internally if not provided.
#' @param m_hat pre-computed estimates on E[Y|X] corresponding to the input x. rboost will compute it internally if not provided.
#' @param ntrees_max the maximum number of trees to grow for xgboost
#' @param num_search_rounds the number of random sampling of hyperparameter combinations for cross validating on xgboost trees
#' @param print_every_n the number of iterations (in each iteration, a tree is grown) by which the code prints out information
#' @param early_stopping_rounds the number of rounds the test error stops decreasing by which the cross validation in finding the optimal number of trees stops
#' @param nthread the number of threads to use. The default is NULL, which uses all available threads
#' @param verbose boolean; whether to print statistic
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rboost_fit = rboost(x, w, y)
#' rboost_est = predict(rboost_fit, x)
#' }
#'
#' @export
rboost= function(x, w, y,
                 k_folds = NULL,
                 p_hat = NULL,
                 m_hat = NULL,
                 ntrees_max = 1000,
                 num_search_rounds = 10,
                 print_every_n = 100,
                 early_stopping_rounds = 10,
                 nthread = NULL,
                 verbose = FALSE){


  input = sanitize_input(x,w,y)
  x = input$x
  w = input$w
  y = input$y
  nobs = nrow(x)
  pobs = ncol(x)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(y)/4)))
  }

  if (is.null(m_hat)){
    y_fit = cvboost(x,
                    y,
                    objective = "reg:squarederror",
                    k_folds = k_folds,
                    ntrees_max = ntrees_max,
                    num_search_rounds = num_search_rounds,
                    print_every_n = print_every_n,
                    early_stopping_rounds = early_stopping_rounds,
                    nthread = nthread,
                    verbose = verbose)

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
                    verbose = verbose)
    p_hat = predict(w_fit)
  }
  else{
    w_fit = NULL
  }

  y_tilde = y - m_hat
  w_tilde = w - p_hat
  pseudo_outcome = y_tilde/w_tilde

  weights = w_tilde^2

  tau_fit = cvboost(x,
                    pseudo_outcome,
                    objective = "reg:squarederror",
                    weights = weights,
                    k_folds = k_folds,
                    ntrees_max = ntrees_max,
                    num_search_rounds = num_search_rounds,
                    print_every_n = print_every_n,
                    early_stopping_rounds = early_stopping_rounds,
                    nthread = nthread,
                    verbose = verbose)

  ret = list(tau_fit = tau_fit,
             pseudo_outcome = pseudo_outcome,
             weights = weights,
             w_fit = w_fit,
             y_fit = y_fit,
             p_hat = p_hat,
             m_hat = m_hat)
  class(ret) = "rboost"
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
#' rboost_fit = rboost(x, w, y)
#' rboost_est = predict(rboost_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.rboost<- function(object,
                          newx = NULL,
                          tau_only = T,
                          ...) {
  if (!is.null(newx)){
    newx = sanitize_x(newx)
  }
  if (tau_only) {
    predict(object$tau_fit, newx = newx)
  } else {
    tau = predict(object$tau_fit, newx = newx)
    e = predict(object$w_fit, newx = newx)
    m = predict(object$y_fit, newx = newx)
    mu1 = m + (1-e) * tau
    mu0 = m - e * tau
    return(list(tau=tau, e=e, m=m, mu1 = mu1, mu0 = mu0))
  }
}
