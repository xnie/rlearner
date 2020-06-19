#' @include utils.R
#'
#' @title S-learner, implemented via kernel ridge regression with a Gaussian kernel
#'
#' @description  S-learner, as proposed by Imai and Ratkovic (2013), implemented via
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds number of folds for cross validation
#' @param b_range the range of Gaussian kernel bandwidths for cross validation
#' @param lambda_range the range of ridge regression penalty factor for cross validation
#'
#' @return an skern object
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' skern_fit = skern(x, w, y)
#' skern_est = predict(skern_fit, x)
#' }
#' @export
skern = function(x, w, y,
                 k_folds = NULL,
                 b_range = 10^(seq(-3,3,0.5)),
                 lambda_range = 10^(seq(-3,3,0.5))){

  c(x, w, y) %<-% sanitize_input(x,w,y)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(w)/4)))
  }

  w = as.numeric(w)
  model_cv= cv_klrs(cbind(x, (w-0.5)*x, (w-0.5)), y, k_folds=k_folds, b_range=b_range, lambda_range=lambda_range)
  mu_1 = predict(model_cv$model, cbind(x, (1-0.5)*x, 1-0.5))$fit
  mu_0 = predict(model_cv$model, cbind(x, (0-0.5)*x, 0-0.5))$fit
  tau_hat = mu_1 - mu_0

  ret = list(tau_fit = model_cv,
             mu_1 = mu_1,
             mu_0 = mu_0,
             tau_hat = tau_hat)
  class(ret) <- "skern"
  ret
}


#' predict for skern
#'
#' get estimated tau(x) using the trained skern model
#'
#' @param object an skern object
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
#' skern_fit = skern(x, w, y)
#' skern_est = predict(skern_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.skern <- function(object,
                          newx = NULL,
                          ...) {
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    mu0_hat = predict(object$tau_fit$model, cbind(newx, (0 - 0.5) * newx, (0 - 0.5)))$fit
    mu1_hat = predict(object$tau_fit$model, cbind(newx, (1 - 0.5) * newx, (1 - 0.5)))$fit
    tau_hat = mu1_hat - mu0_hat
  }
  else {
    tau_hat = object$fit
  }
  return(tau_hat)
}
