#' X-learner, as proposed by Künzel, Sekhon, Bickel, and Yu 2017, implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds_mu1 number of folds for learning E[Y|X,W=1]
#' @param k_folds_mu0 number of folds for learning E[Y|X,W=0]
#' @param k_folds_p number of folds for learning E[W|X]
#' @param lambda.choice how to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param mu1_hat pre-computed estimates on E[Y|X,W=1] corresponding to the input X. xlasso will compute it internally if not provided.
#' @param mu0_hat pre-computed estimates on E[Y|X,W=0] corresponding to the input X. xlasso will compute it internally if not provided.
#' @param p_hat pre-computed estimates on E[W|X] corresponding to the input X. xlasso will compute it internally if not provided
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' xlasso.fit = xlasso(x, w, y)
#' xlasso.est = predict(xlasso.fit, x)
#' }
#'
#'
#' @export
xlasso = function(X, W, Y,
                  alpha=1,
                  k_folds_mu1=NULL,
                  k_folds_mu0=NULL,
                  k_folds_p=NULL,
                  lambda.choice=c("lambda.min", "lambda.1se"),
                  mu1_hat=NULL,
                  mu0_hat=NULL,
                  p_hat=NULL){

  lambda.choice = match.arg(lambda.choice)

  X.1 = X[which(W==1),]
  X.0 = X[which(W==0),]

  Y.1 = Y[which(W==1)]
  Y.0 = Y[which(W==0)]

  nobs.1 = nrow(X.1)
  nobs.0 = nrow(X.0)

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(k_folds_mu1)) {
    k_folds_mu1 = floor(max(3, min(10,nobs.1/4)))
  }

  if (is.null(k_folds_mu0)) {
    k_folds_mu0 = floor(max(3, min(10,nobs.0/4)))
  }

  if (is.null(k_folds_p)) {
    k_folds_p = floor(max(3, min(10,nobs/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid.1 = sample(rep(seq(k_folds_mu1), length = nobs.1))
  foldid.0 = sample(rep(seq(k_folds_mu0), length = nobs.0))
  foldid.W = sample(rep(seq(k_folds_p), length = nobs))

  if (is.null(mu1_hat)){
    t.1.fit = glmnet::cv.glmnet(X.1, Y.1, foldid = foldid.1, alpha = alpha)
    mu1_hat = predict(t.1.fit, newx=X, s=lambda.choice)
  }

  if (is.null(mu0_hat)){
    t.0.fit = glmnet::cv.glmnet(X.0, Y.0, foldid = foldid.0, alpha = alpha)
    mu0_hat = predict(t.0.fit, newx=X, s=lambda.choice)
  }

  D.1 = Y.1 - mu0_hat[W==1]
  D.0 = mu1_hat[W==0] - Y.0

  x.1.fit = glmnet::cv.glmnet(X.1, D.1, foldid = foldid.1, alpha = alpha)
  x.0.fit = glmnet::cv.glmnet(X.0, D.0, foldid = foldid.0, alpha = alpha)

  tau.1.pred = predict(x.1.fit, newx=X, s=lambda.choice)
  tau.0.pred = predict(x.0.fit, newx=X, s=lambda.choice)

  if (is.null(p_hat)){
    w.fit = glmnet::cv.glmnet(X, W, foldid = foldid.W, keep = TRUE, family = "binomial", type.measure = "deviance", alpha = alpha)
    p_hat = w.fit$fit.preval[,!is.na(colSums(w.fit$fit.preval))][, w.fit$lambda == w.fit$lambda.min]
  }
  else{
    w.fit = NULL
  }

  tau.hat = tau.1.pred * (1-p_hat) + tau.0.pred * p_hat

  ret = list(t.1.fit = t.1.fit,
             t.0.fit = t.0.fit,
             x.1.fit = x.1.fit,
             x.0.fit = x.0.fit,
             w.fit = w.fit,
             mu1_hat = mu1_hat,
             mu0_hat = mu0_hat,
             tau.1.pred = tau.1.pred,
             tau.0.pred = tau.0.pred,
             p_hat = p_hat,
             tau.hat = tau.hat)
  class(ret) <- "xlasso"
  ret

}

#' predict for xlasso
#'
#' get estimated tau(x) using the trained xlasso model
#'
#' @param object a xlasso object
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
#' xlasso.fit = xlasso(x, w, y)
#' xlasso.est = predict(xlasso.fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.xlasso <- function(object,
                           newx=NULL,
                           s=c("lambda.min", "lambda.1se"),
                           ...) {
  s = match.arg(s)
  if (!is.null(newx)) {
    tau.1.pred = predict(object$x.1.fit, newx=newx, s=s)
    tau.0.pred = predict(object$x.0.fit, newx=newx, s=s)
    p_hat = predict(object$w.fit, newx=newx, s=s, type="response")
    tau.hat = tau.1.pred * (1-p_hat) + tau.0.pred * p_hat
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
