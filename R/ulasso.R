#' @include utils.R
#'
#' @title U-learner implemented via glmnet (lasso)
#'
#' @description U-learner as proposed by Kunzel, Sekhon, Bickel, and Yu (2017), implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds number of folds for cross-fitting
#' @param lambda_y user-supplied lambda sequence for cross validation in learning E[y|x]
#' @param lambda_w user-supplied lambda sequence for cross validation in learning E[w|x]
#' @param lambda_tau user-supplied lambda sequence for cross validation in learning the treatment effect E[y(1) - y(0) | x]
#' @param lambda_choice how to cross-validate for the treatment effect tau; choose from "lambda.1se" or "lambda.min"
#' @param p_hat user-supplied estimate for E[W|X]
#' @param m_hat user-supplied estimte for E[Y|X]
#' @param cutoff the threshold to cutoff propensity estimate
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' ulasso_fit = ulasso(x, w, y)
#' ulasso_est = predict(ulasso_fit, x)
#' }
#' @export
ulasso = function(x, w, y,
                  alpha = 1,
                  k_folds = NULL,
                  lambda_y = NULL,
                  lambda_w = NULL,
                  lambda_tau = NULL,
                  lambda_choice=c("lambda.1se", "lambda.min"),
                  p_hat = NULL,
                  m_hat = NULL,
                  cutoff = 0.05){

  input = sanitize_input(x,w,y)
  x = input$x
  w = input$w
  y = input$y
  if (!is.logical(w)) {
    stop("w should be a logical vector")
  }

  lambda_choice = match.arg(lambda_choice)

  nobs = nrow(x)
  pobs = ncol(x)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(w)/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid = sample(rep(seq(k_folds), length = length(w)))

  if (is.null(m_hat)){
    y_fit = glmnet::cv.glmnet(x, y, foldid = foldid, keep = TRUE, lambda = lambda_y, alpha = alpha)
    m_hat = y_fit$fit.preval[,!is.na(colSums(y_fit$fit.preval))][, y_fit$lambda == y_fit$lambda.min]
  }
  else {
    y_fit = NULL
  }

    if (is.null(p_hat)){

      w_fit = glmnet::cv.glmnet(x, w,
                               foldid = foldid,
                               family="binomial",
                               type.measure="deviance",
                               keep = TRUE,
                               lambda = lambda_w,
                               alpha = alpha)

      w_lambda_min = w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
      theta_hat = w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
      p_hat = 1/(1 + exp(-theta_hat))
    }
    else{
      w_fit = NULL
    }

  p_hat_thresh = pmax(cutoff, pmin(1 - cutoff, p_hat))

  y_tilde = y - m_hat
  w_tilde = w - p_hat_thresh

  u = y_tilde / w_tilde

  tau_fit = glmnet::cv.glmnet(x, u, nfolds = k_folds, alpha = alpha, lambda = lambda_tau)
  tau_hat = predict(tau_fit, newx = x, s = lambda_choice)

  ret = list(tau_fit = tau_fit,
             w_fit = w_fit,
             y_fit = y_fit,
             p_hat = p_hat,
             m_hat = m_hat,
             tau_hat = tau_hat)
  class(ret) <- "ulasso"
  ret

}

#' predict for ulasso
#'
#' get estimated tau(x) using the trained ulasso model
#'
#' @param object a ulasso object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param ... additional arguments (currently not used)
#' @param s choose from "lambda.min" or "lambda.1se" for prediction
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' ulasso_fit = ulasso(x, w, y)
#' ulasso_est = predict(ulasso_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.ulasso <- function(object,
                           newx = NULL,
                           s = c("lambda.1se", "lambda.min"),
                           ...) {
  s = match.arg(s)
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    tau_hat = predict(object$tau_fit, newx = newx, s = s)
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
