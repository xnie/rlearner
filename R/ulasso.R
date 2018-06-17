#' U-learner, as proposed by Künzel, Sekhon, Bickel, and Yu 2017, implemented via glmnet (lasso)
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the treatment variable (0 or 1)
#' @param alpha tuning parameter for the elastic net
#' @param nfolds number of folds for cross-fitting
#' @param lambda.choice how to cross-validate; choose from "lambda.1se" or "lambda.mse"
#' @param w.hat user-supplied estimate for E[W|X]
#' @param y.hat user-supplied estimte for E[Y|X]
#' @param cutoff the threshold to cutoff propensity estimate
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#'
#' ulasso.fit = ulasso(X, Y, W)
#' ulasso.est = predict(ulasso.fit, X)
#' }
#' @export
ulasso = function(X, Y, W,
                  alpha = 1,
                  nfolds=NULL,
                  lambda.choice=c("lambda.1se", "lambda.mse"),
                  w.hat = NULL,
                  y.hat = NULL,
                  cutoff=0.05){

  lambda.choice = match.arg(lambda.choice)

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,length(W)/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid = sample(rep(seq(nfolds), length = length(W)))

  if (is.null(y.hat)){
    y.fit = glmnet::cv.glmnet(X, Y, foldid = foldid, keep = TRUE, alpha = alpha)
    y.hat = y.fit$fit.preval[,!is.na(colSums(y.fit$fit.preval))][, y.fit$lambda == y.fit$lambda.min]
  }
  else {
    y.fit = NULL
  }

  if (is.null(w.hat)){
    w.fit = glmnet::cv.glmnet(X, W, foldid = foldid, keep = TRUE, family = "binomial", type.measure = "deviance", alpha = alpha)
    w.hat = w.fit$fit.preval[,!is.na(colSums(w.fit$fit.preval))][, w.fit$lambda == w.fit$lambda.min]
  }
  else{
    w.fit = NULL
  }

  w.hat.thresh = pmax(cutoff, pmin(1 - cutoff, w.hat))

  Y.tilde = Y - y.hat
  W.tilde = W - w.hat.thresh

  U = Y.tilde / W.tilde

  tau.fit = glmnet::cv.glmnet(X, U, nfolds= 10, alpha = alpha)
  tau.hat = predict(tau.fit, newx = X, s=lambda.choice)

  ret = list(tau.fit = tau.fit,
             w.fit = w.fit,
             y.fit = y.fit,
             w.hat = w.hat,
             y.hat = y.hat,
             tau.hat = tau.hat)
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
#' X = matrix(rnorm(n*p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[,1], 0) * W + X[,2] + pmin(X[,3], 0) + rnorm(n)
#'
#' ulasso.fit = ulasso(X, Y, W)
#' ulasso.est = predict(ulasso.fit, X)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.ulasso <- function(object,
                           newx=NULL,
                           s=c("lambda.1se", "lambda.min"),
                           ...) {
  s = match.arg(s)
  if (!is.null(newx)) {
    tau.hat = predict(object$tau.fit, newx=newx, s=s)
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
