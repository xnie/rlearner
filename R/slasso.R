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
#' @param lambda user-supplied lambda sequence for cross validation
#' @param lambda_choice how to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty_factor user-supplied penalty factor, must be of length the same as number of features in x
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
                  foldid = NULL,
                  lambda = NULL,
                  lambda_choice = c("lambda.min", "lambda.1se"),
                  penalty_factor = NULL){

  input = sanitize_input(x,w,y)
  x = input$x
  w = input$w
  y = input$y

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

  x_scl_tilde = cbind(as.numeric(2 * w - 1) * cbind(1, x_scl), x_scl)
  x_scl_pred = cbind(1, x_scl, 0 * x_scl)

  if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
    if (!is.null(penalty_factor) && length(penalty_factor) != 2 * pobs + 1) {
      warning("penalty_factor supplied is not 1 plus 2 times the number of columns in x. Using all ones instead.")
    }
    penalty_factor = c(0, rep(1, 2 * pobs))
  }

  s_fit = glmnet::cv.glmnet(x_scl_tilde, y, foldid = foldid, lambda = lambda,
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
