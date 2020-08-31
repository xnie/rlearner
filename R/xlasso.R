#' @include utils.R
#'
#' @title X-learner implemented via glmnet (lasso)
#'
#' @description X-learner as proposed by Kunzel, Sekhon, Bickel, and Yu (2017), implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds_mu1 number of folds for learning E[Y|X,W=1]
#' @param k_folds_mu0 number of folds for learning E[Y|X,W=0]
#' @param k_folds_p number of folds for learning E[W|X]
#' @param lambda_t user-supplied lambda sequence for cross validation in learning E[y|x,w=0] and E[y|x,w=1]
#' @param lambda_x user-supplied lambda sequence for cross validation in learning E[d1|x] and E[d0|x] where d1 = y1 - E[y|x,w=0] and d0 = y0 - E[y|x,w=1]
#' @param lambda_w user-supplied lambda sequence for cross validation in learning E[w|x]
#' @param lambda_choice how to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param mu1_hat pre-computed estimates on E[Y|X,W=1] corresponding to the input x. xlasso will compute it internally if not provided.
#' @param mu0_hat pre-computed estimates on E[Y|X,W=0] corresponding to the input x. xlasso will compute it internally if not provided.
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' xlasso_fit = xlasso(x, w, y)
#' xlasso_est = predict(xlasso_fit, x)
#' }
#'
#'
#' @export
xlasso = function(x, w, y,
                  alpha = 1,
                  k_folds_mu1 = NULL,
                  k_folds_mu0 = NULL,
                  k_folds_p = NULL,
                  lambda_t = NULL,
                  lambda_x = NULL,
                  lambda_w = NULL,
                  lambda_choice = c("lambda.min", "lambda.1se"),
                  mu1_hat = NULL,
                  mu0_hat = NULL){

  c(x, w, y) %<-% sanitize_input(x,w,y)
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

  # fold ID for cross-validation; balance treatment assignments
  foldid_1 = sample(rep(seq(k_folds_mu1), length = nobs_1))
  foldid_0 = sample(rep(seq(k_folds_mu0), length = nobs_0))
  foldid_w = sample(rep(seq(k_folds_p), length = nobs))

  if (is.null(mu1_hat)){
    t_1_fit = glmnet::cv.glmnet(x_1, y_1, foldid = foldid_1, lambda = lambda_t, alpha = alpha)
    mu1_hat = predict(t_1_fit, newx = x, s = lambda_choice)
  }

  if (is.null(mu0_hat)){
    t_0_fit = glmnet::cv.glmnet(x_0, y_0, foldid = foldid_0, lambda = lambda_t, alpha = alpha)
    mu0_hat = predict(t_0_fit, newx = x, s = lambda_choice)
  }

  d_1 = y_1 - mu0_hat[w==1]
  d_0 = mu1_hat[w==0] - y_0

  x_1_fit = glmnet::cv.glmnet(x_1, d_1, foldid = foldid_1, lambda = lambda_x, alpha = alpha)
  x_0_fit = glmnet::cv.glmnet(x_0, d_0, foldid = foldid_0, lambda = lambda_x, alpha = alpha)

  tau_1_pred = predict(x_1_fit, newx = x, s = lambda_choice)
  tau_0_pred = predict(x_0_fit, newx = x, s = lambda_choice)

    w_fit = glmnet::cv.glmnet(x, w,
                             foldid = foldid,
                             family="binomial",
                             type.measure="deviance",
                             lambda = lambda_w,
                             keep = TRUE,
                             alpha = alpha,
                             penalty.factor = penalty_factor_nuisance)

    w_lambda_min = w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
    theta_hat = w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
    p_hat = 1/(1 + exp(-theta_hat))

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
  class(ret) <- "xlasso"
  ret

}

#' predict for xlasso
#'
#' get estimated tau(x) using the trained xlasso model
#'
#' @param object a xlasso object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param s choose from "lambda.min" or "lambda.1se
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
#' xlasso_fit = xlasso(x, w, y)
#' xlasso_est = predict(xlasso_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.xlasso <- function(object,
                           newx = NULL,
                           s = c("lambda.min", "lambda.1se"),
                           ...) {
  s = match.arg(s)
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    tau_1_pred = predict(object$x_1_fit, newx = newx, s = s)
    tau_0_pred = predict(object$x_0_fit, newx = newx, s = s)
    p_hat = predict(object$w_fit, newx=newx, s = s, type = "response")
    tau_hat = tau_1_pred * (1 - p_hat) + tau_0_pred * p_hat
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
