#' @include utils.R
#'
#' @title T-learner, implemented via glmnet (lasso)
#'
#' @description T-learner learns the treated and control expected outcome respectively by fitting two separate models.
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds_mu1 number of folds for cross validation for the treated
#' @param k_folds_mu0 number of folds for cross validation for the control
#' @param lambda user-supplied lambda sequence for cross validation
#' @param lambda_choice how to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty_factor user-supplied penalty factor, must be of length the same as number of features in x
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0_5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' tlasso_fit = tlasso(x, w, y)
#' tlasso_est = predict(tlasso_fit, x)
#' }
#' @export
tlasso = function(x, w, y,
                  alpha = 1,
                  k_folds_mu1 = NULL,
                  k_folds_mu0 = NULL,
                  lambda = NULL,
                  lambda_choice = c("lambda.min", "lambda.1se"),
                  penalty_factor= NULL) {


  input = sanitize_input(x,w,y)
  x = input$x
  w = input$w
  y = input$y
  if (!is.logical(w)) {
    stop("w should be a logical vector")
  }


  lambda_choice = match.arg(lambda_choice)

  x_1 = x[which(w==1),]
  x_0 = x[which(w==0),]

  y_1 = y[which(w==1)]
  y_0 = y[which(w==0)]

  nobs_1 = nrow(x_1)
  nobs_0 = nrow(x_0)

  pobs = ncol(x)

  if (is.null(k_folds_mu1)) {
    k_folds_mu1 = floor(max(3, min(10, nobs_1/4)))
  }

  if (is.null(k_folds_mu0)) {
    k_folds_mu0 = floor(max(3, min(10, nobs_0/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid_1 = sample(rep(seq(k_folds_mu1), length = nobs_1))
  foldid_0 = sample(rep(seq(k_folds_mu0), length = nobs_0))

  if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
    penalty_factor = rep(1, pobs)
    if (!is.null(penalty_factor) && length(penalty_factor) != pobs) {
      warning("penalty_factor supplied is not of the same length as the number of columns in x. Using all ones instead.")
    }
  }

  t_1_fit = glmnet::cv.glmnet(x_1, y_1, foldid = foldid_1, alpha = alpha, lambda = lambda, penalty.factor=penalty_factor)
  t_0_fit = glmnet::cv.glmnet(x_0, y_0, foldid = foldid_0, alpha = alpha, lambda = lambda, penalty.factor=penalty_factor)

  y_1_pred = predict(t_1_fit, newx = x, s = lambda_choice)
  y_0_pred = predict(t_0_fit, newx = x, s = lambda_choice)

  tau_hat = y_1_pred - y_0_pred

  ret = list(t_1_fit = t_1_fit,
             t_0_fit = t_0_fit,
             y_1_pred = y_1_pred,
             y_0_pred = y_0_pred,
             tau_hat = tau_hat)
  class(ret) <- "tlasso"
  ret
}

#' predict for tlasso
#'
#' get estimated tau(x) using the trained tlasso model
#'
#' @param object a tlasso object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param s choose from "lambda.min" or "lambda.1se" for prediction
#' @param ... additional arguments (currently not used)
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0_5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' tlasso_fit = tlasso(x, w, y)
#' tlasso_est = predict(tlasso_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.tlasso <- function(object,
                           newx = NULL,
                           s = c("lambda.min", "lambda.1se"),
                           ...) {
  s = match.arg(s)

  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    y_1_pred = predict(object$t_1_fit, newx = newx, s = s)
    y_0_pred = predict(object$t_0_fit, newx = newx, s = s)
    tau_hat = y_1_pred - y_0_pred
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
