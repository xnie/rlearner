#' @include utils.R
#'
#' @title S-learner, implemented via glmnet (lasso)
#'
#' @description  S-learner, as proposed by Imai and Ratkovic (2013), implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds number of folds for cross validation
#' @param lambda_choice how to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty_search whether to perform fine grainted penalty_factor search (logical)
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' slasso_fit = slasso(x, w, y)
#' slasso_est = predict(slasso_fit, x)
#' }
#' @export
slasso = function(x, w, y,
                  alpha = 1,
                  k_folds = NULL,
                  lambda_choice = c("lambda.min", "lambda.1se"),
                  penalty_search = FALSE) {

  c(x, w, y) %<-% sanitize_input(x,w,y)

  standardization = caret::preProcess(x, method=c("center", "scale")) # get the standardization params
  x_scl = predict(standardization, x)							 # standardize the input
  x_scl = x_scl[,!is.na(colSums(x_scl)), drop = FALSE]

  lambda_choice = match.arg(lambda_choice)

  nobs = nrow(x_scl)
  pobs = ncol(x_scl)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10, nobs/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid = sample(rep(seq(k_folds), length = nobs))

  x_scl_tilde = cbind(as.numeric(2 * w - 1) * cbind(1, x_scl), x_scl)
  x_scl_pred = cbind(1, x_scl, 0 * x_scl)
  if (penalty_search) {
    search_range = 5
    cvm_min = Inf
    last_best = NULL
    for (l in 0:2){

      updated = FALSE

      for (i in 1:search_range){

        if (l==0){
          penalty_factor = c(0, rep(10^(i- ceiling(search_range / 2.0)), pobs), rep(1, pobs))
        }
        else{
          penalty_factor = c(0, rep(10^(last_best - 10^(-l + 1) + 20.0 / search_range * 10^(-l) * i), pobs), rep(1, pobs))
        }
        s_fit <- glmnet::cv.glmnet(x = x_scl_tilde, y = y, foldid = foldid, penalty.factor = penalty_factor, standardize = FALSE, alpha = alpha)
        s_fit_cvm = s_fit$cvm[s_fit$lambda == s_fit$lambda.min]
        s_beta = as.vector(t(coef(s_fit, s=lambda_choice)[-1]))
        if (s_fit_cvm < cvm_min){
          cvm_min = s_fit_cvm
          s_fit_best = s_fit
          if (l==0){
            best_i = i - ceiling(search_range/2.0)
          }
          else{
            best_i = i
          }
          best_penalty_factor = penalty_factor
          updated = TRUE
        }
      }

      if (l==0){
        last_best = best_i
      }
      else{
        if (updated){
          last_best = last_best - 10^(-l + 1) + 20 / search_range * 10^(-l) * best_i
        }
      }
    }
    penalty_factor = best_penalty_factor
  }
  else{
    penalty_factor = c(0, rep(1, 2 * pobs))
  }

  s_fit = glmnet::cv.glmnet(x_scl_tilde, y, foldid = foldid,
                            penalty.factor = penalty_factor,
                            standardize = FALSE, alpha = alpha)

  s_beta = as.vector(t(coef(s_fit, s = lambda_choice)[-1]))

  tau_hat = 2 * x_scl_pred %*% s_beta

  ret = list(s_fit = s_fit,
             s_beta = s_beta,
             tau_hat = tau_hat,
             standardization = standardization)

  class(ret) <- "slasso"
  ret
}

#' predict for slasso
#'
#' get estimated tau(x) using the trained slasso model
#'
#' @param object a slasso object
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
#' slasso_fit = slasso(x, w, y)
#' slasso_est = predict(slasso_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.slasso <- function(object,
                           newx = NULL,
                           ...) {
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    newx_scl = predict(object$standardization, newx) # standardize the new data using the same standardization as with the training data
    newx_scl = newx_scl[,!is.na(colSums(newx_scl)), drop = FALSE]
    newx_scl_pred = cbind(1, newx_scl, 0 * newx_scl)
    tau_hat = 2 * newx_scl_pred %*% object$s_beta
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
