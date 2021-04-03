#' @include utils.R
#'
#' @title X-learner implemented via kernel ridge regression (with a Gaussian kernel)
#'
#' @description X-learner as proposed by Kunzel, Sekhon, Bickel, and Yu (2017), implemented via kernel ridge regression (with a Gaussian kernel)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds number of folds for cross-fitting
#' @param b_range the range of Gaussian kernel bandwidths for cross validation
#' @param lambda_range the range of ridge regression penalty factor for cross validation
#' @param mu1_hat pre-computed estimates on E[Y|X,W=1] corresponding to the input x. xkern will compute it internally if not provided.
#' @param mu0_hat pre-computed estimates on E[Y|X,W=0] corresponding to the input x. xkern will compute it internally if not provided.
#' @param p_hat user-supplied estimate for E[W|X]. xkern will compute it internally if not provided.
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' xkern_fit = xkern(x, w, y)
#' xkern_est = predict(xkern_fit, x)
#' }
#'
#' @export
xkern= function(x, w, y,
                k_folds=NULL,
                b_range =10^(seq(-3,3,0.5)),
                lambda_range = 10^(seq(-3,3,0.5)),
                mu1_hat=NULL,
                mu0_hat=NULL,
                p_hat=NULL){


  input = sanitize_input(x,w,y)
  x = input$x
  w = input$w
  y = input$y
  w = as.numeric(w)
  x_1 = x[which(w==1),]
  x_0 = x[which(w==0),]

  y_1 = y[which(w==1)]
  y_0 = y[which(w==0)]

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(w)/4)))
  }

  if (is.null(mu1_hat)){
    t_1_fit = cv_klrs(x_1, y_1, k_folds=k_folds, b_range=b_range, lambda_range=lambda_range)
    mu1_hat = predict(t_1_fit$model, x)$fit
  } else{
    t_1_fit = NULL
  }
  if (is.null(mu0_hat)){
    t_0_fit = cv_klrs(x_0, y_0, k_folds=k_folds, b_range=b_range, lambda_range=lambda_range)
    mu0_hat = predict(t_0_fit$model, x)$fit
  } else {
    t_0_fit = NULL
  }

  d_1 = y_1 - mu0_hat[w==1]
  d_0 = mu1_hat[w==0] - y_0


  x_1_fit = cv_klrs(x_1, d_1, k_folds=k_folds, b_range=b_range, lambda_range=lambda_range)
  x_0_fit = cv_klrs(x_0, d_0, k_folds=k_folds, b_range=b_range, lambda_range=lambda_range)

  tau_1_pred = predict(x_1_fit$model, x)$fit
  tau_0_pred = predict(x_0_fit$model, x)$fit

  if (is.null(p_hat)) {
   p_hat_model = cv_klrs(x, w, weights=NULL, k_folds=k_folds, b_range=b_range,lambda_range=lambda_range)
   p_hat = p_hat_model$fit
  } else {
   p_hat_model = NULL
  }

  tau_hat = tau_1_pred * (1 - p_hat) + tau_0_pred * p_hat

  ret = list(t_1_fit = t_1_fit,
             t_0_fit = t_0_fit,
             x_1_fit = x_1_fit,
             x_0_fit = x_0_fit,
             p_hat_model = p_hat_model,
             mu1_hat = mu1_hat,
             mu0_hat = mu0_hat,
             tau_1_pred = tau_1_pred,
             tau_0_pred = tau_0_pred,
             p_hat = p_hat,
             tau_hat = tau_hat)
  class(ret) <- "xkern"
  ret
}

#' predict for xkern
#'
#' get estimated tau(x) using the trained xkern model
#'
#' @param object an xkern object
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
#' xkern_fit = xkern(x, w, y)
#' xkern_est = predict(xkern_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.xkern<- function(object,
                           newx = NULL,
                           ...) {
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    tau_1_pred = predict(object$x_1_fit$model, newx)$fit
    tau_0_pred = predict(object$x_0_fit$model, newx)$fit
    p_hat = predict(object$p_hat_model$model, newx)$fit
    tau_hat = tau_1_pred * (1 - p_hat) + tau_0_pred * p_hat
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
