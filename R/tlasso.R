#' T-learner, implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds_mu1 number of folds for cross validation for the treated
#' @param k_folds_mu0 number of folds for cross validation for the control
#' @param lambda.choice how to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' tlasso.fit = tlasso(x, w, y)
#' tlasso.est = predict(tlasso.fit, x)
#' }
#' @export
tlasso = function(X, W, Y,
                  alpha = 1,
                  k_folds_mu1=NULL,
                  k_folds_mu0=NULL,
                  lambda.choice=c("lambda.min", "lambda.1se")) {

  lambda.choice = match.arg(lambda.choice)

  X.1 = X[which(W==1),]
  X.0 = X[which(W==0),]

  Y.1 = Y[which(W==1)]
  Y.0 = Y[which(W==0)]

  nobs.1 = nrow(X.1)
  nobs.0 = nrow(X.0)

  pobs = ncol(X)

  if (is.null(k_folds_mu1)) {
    k_folds_mu1 = floor(max(3, min(10,nobs.1/4)))
  }

  if (is.null(k_folds_mu0)) {
    k_folds_mu0 = floor(max(3, min(10,nobs.0/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid.1 = sample(rep(seq(k_folds_mu1), length = nobs.1))
  foldid.0 = sample(rep(seq(k_folds_mu0), length = nobs.0))

  t.1.fit = glmnet::cv.glmnet(X.1, Y.1, foldid = foldid.1, alpha = alpha)
  t.0.fit = glmnet::cv.glmnet(X.0, Y.0, foldid = foldid.0, alpha = alpha)

  y.1.pred = predict(t.1.fit, newx=X, s=lambda.choice)
  y.0.pred = predict(t.0.fit, newx=X, s=lambda.choice)

  tau.hat = y.1.pred - y.0.pred

  ret = list(t.1.fit = t.1.fit,
             t.0.fit = t.0.fit,
             y.1.pred = y.1.pred,
             y.0.pred = y.0.pred,
             tau.hat = tau.hat)
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
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' tlasso.fit = tlasso(x, w, y)
#' tlasso.est = predict(tlasso.fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.tlasso <- function(object,
                           newx=NULL,
                           s=c("lambda.min", "lambda.1se"),
                           ...) {
  s = match.arg(s)

  if (!is.null(newx)) {
    y.1.pred = predict(object$t.1.fit, newx=newx, s=s)
    y.0.pred = predict(object$t.0.fit, newx=newx, s=s)
    tau.hat = y.1.pred - y.0.pred
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
