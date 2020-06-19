#' @include utils.R
#'
#' @title T-learner, implemented via kernel ridge regression with Gaussian kernels
#'
#' @description T-learner learns the treated and control expected outcome respectively by fitting two separate models.
#'
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param b_range the range of Gaussian kernel bandwidths for cross validation
#' @param lambda_range the range of ridge regression penalty factor for cross validation
#' @param k_folds number of folds for cross validation
#'
#' @return an tkern object
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' tkern_fit = tkern(x, w, y)
#' tkern_est = predict(tkern_fit, x)
#' }
#' @export
tkern = function(x, w, y,
                 b_range =10^(seq(-3,3,0.5)),
                 lambda_range = 10^(seq(-3,3,0.5)),
                 k_folds = NULL){

  c(x, w, y) %<-% sanitize_input(x,w,y)

  w = as.numeric(w)

  x_1 = x[which(w==1),]
  x_0 = x[which(w==0),]

  y_1 = y[which(w==1)]
  y_0 = y[which(w==0)]

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(w)/4)))
  }

  t_1_fit = cv_klrs(x_1, y_1, k_folds=k_folds, b_range=b_range, lambda_range=lambda_range)
  t_0_fit = cv_klrs(x_0, y_0, k_folds=k_folds, b_range=b_range, lambda_range=lambda_range)

  y_1_pred = predict(t_1_fit$model, x)$fit
  y_0_pred = predict(t_0_fit$model, x)$fit

  tau_hat = y_1_pred - y_0_pred

  ret = list(t_1_fit = t_1_fit,
             t_0_fit = t_0_fit,
             y_1_pred = y_1_pred,
             y_0_pred = y_0_pred,
             tau_hat = tau_hat)
  class(ret) <- "tkern"
  ret
}


#' predict for tkern
#'
#' get estimated tau(x) using the trained tkern model
#'
#' @param object an tkern object
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
#' tkern_fit = tkern(x, w, y)
#' tkern_est = predict(tkern_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.tkern <- function(object,
                           newx = NULL,
                           ...) {
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    y_1_pred = predict(object$t_1_fit$model, newx)$fit
    y_0_pred = predict(object$t_0_fit$model, newx)$fit
    tau_hat = y_1_pred - y_0_pred
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
