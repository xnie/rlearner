#' X-lasso as proposed by Kunzel et al 2017
#'
#' @param X
#' @param Y
#' @param W
#' @param alpha
#' @param nfolds.1
#' @param nfolds.0
#' @param nfolds.W
#' @param lambda.choice
#'
#' @return
#' @export xlasso
#'
#' @examples
xlasso = function(X, Y, W,
                  alpha = 1,
                  nfolds.1=NULL,
                  nfolds.0=NULL,
                  nfolds.W=NULL,
                  lambda.choice=c("lambda.min", "lambda.1se")) {

  lambda.choice = match.arg(lambda.choice)

  X.1 = X[which(W==1),]
  X.0 = X[which(W==0),]

  Y.1 = Y[which(W==1)]
  Y.0 = Y[which(W==0)]

  nobs.1 = nrow(X.1)
  nobs.0 = nrow(X.0)

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(nfolds.1)) {
    nfolds.1 = floor(max(3, min(10,nobs.1/4)))
  }

  if (is.null(nfolds.0)) {
    nfolds.0 = floor(max(3, min(10,nobs.0/4)))
  }

  if (is.null(nfolds.W)) {
    nfolds.W = floor(max(3, min(10,nobs/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid.1 = sample(rep(seq(nfolds.1), length = nobs.1))
  foldid.0 = sample(rep(seq(nfolds.0), length = nobs.0))
  foldid.W = sample(rep(seq(nfolds.W), length = nobs))

  t.1.fit = glmnet::cv.glmnet(X.1, Y.1, foldid = foldid.1)
  t.0.fit = glmnet::cv.glmnet(X.0, Y.0, foldid = foldid.0)

  y.1.pred = predict(t.1.fit, newx=X, s=lambda.choice)
  y.0.pred = predict(t.0.fit, newx=X, s=lambda.choice)

  D.1 = Y.1 - y.0.pred[W==1]
  D.0 = y.1.pred[W==0] - Y.0

  x.1.fit = glmnet::cv.glmnet(X.1, D.1, foldid = foldid.1)
  x.0.fit = glmnet::cv.glmnet(X.0, D.0, foldid = foldid.0)

  tau.1.pred = predict(x.1.fit, newx=X, s=lambda.choice)
  tau.0.pred = predict(x.0.fit, newx=X, s=lambda.choice)

  w.fit = glmnet::cv.glmnet(X, W, foldid=foldid.W, keep=TRUE, family="binomial", type.measure = "auc", alpha = alpha)
  w.hat = w.fit$fit.preval[,!is.na(colSums(w.fit$fit.preval))][, w.fit$lambda == w.fit$lambda.min]


  tau.hat = tau.1.pred * (1-w.hat) + tau.0.pred * w.hat

  ret = list(t.1.fit = t.1.fit,
             t.0.fit = t.0.fit,
             x.1.fit = x.1.fit,
             x.0.fit = x.0.fit,
             w.fit = w.fit,
             y.1.pred = y.1.pred,
             y.0.pred = y.0.pred,
             tau.1.pred = tau.1.pred,
             tau.0.pred = tau.0.pred,
             w.hat = w.hat,
             tau.hat = tau.hat)
  class(ret) <- "xlasso"
  ret

}

predict.xlasso <- function(object,
                           newx=NULL,
                           s=c("lambda.min", "lambda.1se"),
                           ...) {
  s = match.arg(s)
  if (!is.null(newx)) {
    tau.1.pred = predict(object$x.1.fit, newx=newx, s=s)
    tau.0.pred = predict(object$x.0.fit, newx=newx, s=s)
    w.hat = predict(object$w.fit, newx=newx, s=s, type="response")
    tau.hat = tau.1.pred * (1-w.hat) + tau.0.pred * w.hat
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
