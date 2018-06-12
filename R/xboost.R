#' X-learner, as proposed by Künzel, Sekhon, Bickel, and Yu 2017
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the treatment variable (0 or 1)
#' @param nfolds.1 number of folds for learning E[Y|X,W=1]
#' @param nfolds.0 number of folds for learning E[Y|X,W=0]
#' @param nfolds.W number of folds for learning E[W|X]
#'
#' @export xboost
xboost = function(X, Y, W,
                  nfolds.1=NULL,
                  nfolds.0=NULL,
                  nfolds.W=NULL){

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

  t.1.fit = cvboost(X.1, Y.1, objective="reg:linear", nfolds = nfolds.1)
  t.0.fit = cvboost(X.0, Y.0, objective="reg:linear", nfolds = nfolds.0)

  y.1.pred = predict(t.1.fit, newx=X)
  y.0.pred = predict(t.0.fit, newx=X)

  D.1 = Y.1 - y.0.pred[W==1]
  D.0 = y.1.pred[W==0] - Y.0

  x.1.fit = cvboost(X.1, D.1, objective="reg:linear", nfolds = nfolds.1)
  x.0.fit = cvboost(X.0, D.0, objective="reg:linear", nfolds = nfolds.0)

  tau.1.pred = predict(x.1.fit, newx=X)
  tau.0.pred = predict(x.0.fit, newx=X)

  w.fit = cvboost(X, W, objective="binary:logistic", nfolds = nfolds.W)
  w.hat = predict(w.fit)

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
  class(ret) <- "xboost"
  ret

}

#' Title
#'
#' @param object
#' @param newx
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict.xboost <- function(object,
                           newx=NULL,
                           ...) {
  if (!is.null(newx)) {
    tau.1.pred = predict(object$x.1.fit, newx=newx)
    tau.0.pred = predict(object$x.0.fit, newx=newx)
    w.hat = predict(object$w.fit, newx=newx)
    tau.hat = tau.1.pred * (1-w.hat) + tau.0.pred * w.hat
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
