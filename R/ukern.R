#' @include utils.R
#'
#' @title U-learner implemented via kernel ridge regression with a Gaussian kernel
#'
#' @description U-learner as proposed by Kunzel, Sekhon, Bickel, and Yu (2017), implemented via kernel ridge regression with a Gaussian kernel
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds number of folds for cross-fitting
#' @param p_hat user-supplied estimate for E[W|X]
#' @param m_hat user-supplied estimte for E[Y|X]
#' @param b_range the range of Gaussian kernel bandwidths for cross validation
#' @param lambda_range the range of ridge regression penalty factor for cross validation
#' @param cutoff the threshold to cutoff propensity estimate
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' ukern_fit = ukern(x, w, y)
#' ukern_est = predict(ukern_fit, x)
#' }
#' @export
ukern = function(x, w, y,
                 k_folds = NULL,
                 p_hat = NULL,
                 m_hat = NULL,
                 b_range = 10^(seq(-3,3,0.5)),
                 lambda_range = 10^(seq(-3,3,0.5)),
                 cutoff = 0.05){

  input = sanitize_input(x,w,y)
  x = input$x
  w = input$w
  y = input$y

  w = as.numeric(w)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(w)/4)))
  }

  if (is.null(p_hat)) {
  p_hat_model = cv_klrs(x, w, weights=NULL, k_folds=k_folds, b_range=b_range,lambda_range=lambda_range)
  p_hat = p_hat_model$fit
 } else {
  p_hat_model = NULL
 }
 if (is.null(m_hat)) {
  m_hat_model = cv_klrs(x, y, weights=NULL, k_folds=k_folds, b_range=b_range,lambda_range=lambda_range)
  m_hat = m_hat_model$fit
 } else {
  m_hat_model = NULL
 }

  p_hat = pmax(cutoff, pmin(1 - cutoff, p_hat))
  y_tilde = y - m_hat
  w_tilde = w - p_hat
  pseudo_outcome = y_tilde/w_tilde

  model_tau_cv= cv_klrs(x, pseudo_outcome, k_folds=k_folds, b_range=b_range,lambda_range=lambda_range)# lambda is the same as var from kernlab. 1/b is the same as sigma in kernlab

  ret = list(tau_fit = model_tau_cv,
             p_fit = p_hat_model,
             m_fit = m_hat_model,
             p_hat = p_hat,
             m_hat = m_hat)
  class(ret) <- "ukern"
  ret
}

#' predict for ukern
#'
#' get estimated tau(x) using the trained ukern model
#'
#' @param object a ukern object
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
#' ukern_fit = ukern(x, w, y)
#' ukern_est = predict(ukern_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.ukern <- function(object,
                          newx = NULL,
                          ...) {
  if (!is.null(newx)){
    newx = sanitize_x(newx)
    tau_hat = predict(object$tau_fit$model, newx)$fit
  } else {
    tau_hat = object$tau_fit$fit
  }
  tau_hat
}

