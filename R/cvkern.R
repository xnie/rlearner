#' Kernel ridge regression with Gaussian kernels with cross validation to search for hyper-parameters (implemented with th eKLRS package)
#'
#' @param x the input features
#' @param y the observed response (real valued)
#' @param weights weights for input if doing weighted regression/classification. If set to NULL, no weights are used
#' @param k_folds number of folds used in cross validation
#' @param b_range the range of Gaussian kernel bandwidths for cross validation
#' @param lambda_range the range of ridge regression penalty factor for cross validation
#'
#' @return a list that includes the best b, lambda, best predictions, best model fit
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' y = pmax(x[,1], 0) + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' fit = cv_klrs(x, y)
#' }
#'
#' @export
cv_klrs = function(x,
                   y,
                   weights=NULL,
                   k_folds=NULL,
                   b_range = 10^(seq(-3,3,0.5)),
                   lambda_range = 10^(seq(-3,3,0.5))) {
  krls_available = requireNamespace("KRLS2", quietly = TRUE)
  if (!krls_available) {
    stop("KRLS2 needs to be installed for kernel ridge regression.")
  } else{
    if (!packageVersion("KRLS2")>"1.1.0") {
      stop("KRLS2 needs to be of version at least 1.1.1.")
    }
  }
  if (!is.null(weights)) {
    truncate = T
  } else {
    truncate = F
  }

  best_mse = NULL

  foldid = sample(rep(seq(k_folds), length = length(y)))
  for (b in b_range) {
    for (lambda in lambda_range) {
      fit = rep(0, length(y))
      for (f in 1:k_folds) {
        model_f = suppressMessages(KRLS2::krls(X=x[foldid != f,], y=y[foldid != f], w=weights[foldid != f], loss="leastsquares", b=b, lambda=lambda, truncate = truncate))# lambda is the same as var from kernlab. 1/b is the same as sigma in kernlab
        fit[foldid == f] = predict(model_f, x[foldid == f,])$fit
      }
      if (is.null(weights)) {
        mse_fit = mean((y-fit)^2)
      } else {
        mse_fit = mean(weights*(y-fit)^2)
      }
      if (is.null(best_mse)) {
        best_b = b
        best_lambda = lambda
        best_mse = mse_fit
        best_fit = fit
      } else if (best_mse > mse_fit) {
        best_b = b
        best_lambda = lambda
        best_mse = mse_fit
        best_fit = fit
      }
    }
  }
  model = suppressMessages(KRLS2::krls(X=x, y=y, w= weights, loss="leastsquares", b=best_b, lambda=best_lambda, truncate = truncate))# lambda is the same as var from kernlab. 1/b is the same as sigma in kernlab
  return(list(b=best_b,
              lambda = best_lambda,
              fit = best_fit,
              model = model
  ))
}
