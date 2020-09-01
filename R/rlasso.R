#' @include utils.R
#'
#' @title R-learner, implemented via glmnet (lasso)
#'
#' @description  R-learner, as proposed by Nie and Wager (2017), implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds number of folds for cross-fitting
#' @param foldid user-supplied foldid. Must have length equal to length(w). If provided, it overrides the k_folds option.
#' @param lambda_y user-supplied lambda sequence for cross validation in learning E[y|x]
#' @param lambda_w user-supplied lambda sequence for cross validation in learning E[w|x]
#' @param lambda_tau user-supplied lambda sequence for cross validation in learning the treatment effect E[y(1) - y(0) | x]
#' @param lambda_choice how to cross-validate for learning the treatment effect tau; choose from "lambda.min" or "lambda.1se"
#' @param rs whether to use the RS-learner (logical).
#' @param p_hat user-supplied estimate for E[W|X]
#' @param m_hat user-supplied estimte for E[Y|X]
#' @param penalty_factor user-supplied penalty factor, a vector of length the same as the number of covariates in x.
#' @return an rlasso object
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rlasso_fit = rlasso(x, w, y)
#' rlasso_est = predict(rlasso_fit, x)
#' }
#' @export
rlasso = function(x, w, y,
                  alpha = 1,
                  k_folds = NULL,
                  foldid = NULL,
                  lambda_y = NULL,
                  lambda_w = NULL,
                  lambda_tau = NULL,
                  lambda_choice = c("lambda.min","lambda.1se"),
                  rs = FALSE,
                  p_hat = NULL,
                  m_hat = NULL,
                  penalty_factor = NULL){

    c(x, w, y) %<-% sanitize_input(x,w,y)

    standardization = caret::preProcess(x, method=c("center", "scale")) # get the standardization params
    x_scl = predict(standardization, x)							 # standardize the input
    x_scl = x_scl[,!is.na(colSums(x_scl)), drop = FALSE]

    lambda_choice = match.arg(lambda_choice)

    nobs = nrow(x_scl)
    pobs = ncol(x_scl)

    if (is.null(foldid) || length(foldid) != length(w)) {

      if (!is.null(foldid) && length(foldid) != length(w)) {
        warning("supplied foldid does not have the same length ")
      }

      if (is.null(k_folds)) {
          k_folds = floor(max(3, min(10,length(w)/4)))
      }

      # fold ID for cross-validation; balance treatment assignments
      foldid = sample(rep(seq(k_folds), length = length(w)))

    }

    # penalty factor for nuisance and tau estimators
    if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
      if (!is.null(penalty_factor) && length(penalty_factor) != pobs) {
        warning("penalty_factor supplied is not of the same length as the number of columns in x after removing NA columns. Using all ones instead.")
      }
      penalty_factor_nuisance = rep(1, pobs)
      if (rs) {
        penalty_factor_tau = c(0, rep(1, 2 * pobs))
      }
      else {
        penalty_factor_tau = c(0, rep(1, pobs))
      }
    }
    else {
      penalty_factor_nuisance = penalty_factor
      if (rs) {
        penalty_factor_tau = c(0, penalty_factor, penalty_factor)
      }
      else {
        penalty_factor_tau = c(0, penalty_factor)
      }
    }

    if (is.null(m_hat)){
      y_fit = glmnet::cv.glmnet(x, y,
                                foldid = foldid,
                                keep = TRUE,
                                lambda = lambda_y,
                                alpha = alpha,
                                penalty.factor = penalty_factor_nuisance)

      y_lambda_min = y_fit$lambda[which.min(y_fit$cvm[!is.na(colSums(y_fit$fit.preval))])]
      m_hat = y_fit$fit.preval[,!is.na(colSums(y_fit$fit.preval))][, y_fit$lambda[!is.na(colSums(y_fit$fit.preval))] == y_lambda_min]
    }
    else {
      y_fit = NULL
    }

    if (is.null(p_hat)){

    if (is.logical(w)) {
      w_fit = glmnet::cv.glmnet(x, w,
                               foldid = foldid,
                               family="binomial",
                               type.measure="deviance",
                               keep = TRUE,
                               lambda = lambda_w,
                               alpha = alpha,
                               penalty.factor = penalty_factor_nuisance)

      w_lambda_min = w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
      theta_hat = w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
      p_hat = 1/(1 + exp(-theta_hat))
  	} else {
      w_fit = glmnet::cv.glmnet(x, w,
                               foldid = foldid,
                               lambda = lambda_w,
                               keep = TRUE,
                               alpha = alpha,
                               penalty.factor = penalty_factor_nuisance)

      w_lambda_min = w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
      p_hat = w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
  	}
    }
    else{
      w_fit = NULL
    }

    y_tilde = y - m_hat

    if (rs){
      x_scl_tilde = cbind(as.numeric(w - p_hat) * cbind(1, x_scl), x_scl)
      x_scl_pred = cbind(1, x_scl, x_scl * 0)
    }
    else{
      x_scl_tilde = cbind(as.numeric(w - p_hat) * cbind(1, x_scl))
      x_scl_pred = cbind(1, x_scl)
    }

    tau_fit = glmnet::cv.glmnet(x_scl_tilde,
                                y_tilde,
                                foldid = foldid,
                                alpha = alpha,
                                lambda = lambda_tau,
                                penalty.factor = penalty_factor_tau,
                                standardize = FALSE)

    tau_beta = as.vector(t(coef(tau_fit, s = lambda_choice)[-1]))

    tau_hat = x_scl_pred %*% tau_beta

    ret = list(tau_fit = tau_fit,
               tau_beta = tau_beta,
               w_fit = w_fit,
               y_fit = y_fit,
               p_hat = p_hat,
               m_hat = m_hat,
               tau_hat = tau_hat,
               rs = rs,
               standardization = standardization)
    class(ret) <- "rlasso"
    ret
}


#' predict for rlasso
#'
#' get estimated tau(x) using the trained rlasso model
#'
#' @param object an rlasso object
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
#' rlasso_fit = rlasso(x, w, y)
#' rlasso_est = predict(rlasso_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.rlasso <- function(object,
                           newx = NULL,
                           ...) {
  if (!is.null(newx)) {

    newx = sanitize_x(newx)
    newx_scl = predict(object$standardization, newx) # standardize the new data using the same standardization as with the training data
    newx_scl = newx_scl[,!is.na(colSums(newx_scl)), drop = FALSE]

    if (object$rs){
      newx_scl_pred = cbind(1, newx_scl, newx_scl * 0)
    }
    else{
      newx_scl_pred = cbind(1, newx_scl)
    }
    tau_hat = newx_scl_pred %*% object$tau_beta
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
