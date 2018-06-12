#' Title
#'
#' @param X
#' @param Y
#' @param W
#' @param nfolds
#' @param rc
#' @param w.hat
#' @param y.hat
#'
#' @return
#' @export
#'
#' @examples
rboost= function(X, Y, W,
                  nfolds=NULL,
                  rc = FALSE,
                  w.hat = NULL,
                  y.hat = NULL,
                  nthread = NULL){

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,length(W)/4)))
  }

  if (is.null(y.hat)){
    y.fit = cvboost(X, Y, objective="reg:linear", nfolds=nfolds)
    y.hat = predict(y.fit)
  }
  else {
    y.fit = NULL
  }

  if (is.null(w.hat)){
    w.fit = cvboost(X, W, objective="binary:logistic", nfolds=nfolds)
    w.hat = predict(w.fit)
  }
  else{
    w.fit = NULL
  }

  Y.tilde = Y - y.hat
  W.tilde = W - w.hat
  if (rc) {
    tau.const.fit = lm(Y.tilde ~ W.tilde)
    tau.const = coef(tau.const.fit)["W.tilde"]
    Y.tilde.tilde = Y.tilde - W.tilde * tau.const # subtracting out the constant treatment effect
    pseudo.outcome = Y.tilde.tilde/W.tilde
  }
  else{
    pseudo.outcome = Y.tilde/W.tilde
    tau.const = NULL
  }

  weights = W.tilde^2

  tau.fit = cvboost(X, pseudo.outcome, objective="reg:linear", weights=weights, nfolds=nfolds)

  ret = list(tau.fit = tau.fit,
             w.fit = w.fit,
             y.fit = y.fit,
             w.hat = w.hat,
             y.hat = y.hat,
             tau.const = tau.const,
             rc = rc)
  class(ret) <- "rboost"
  ret
}

#' Title
#'
#' @param object
#' @param newx
#' @param ...
#'
#' @return
#' @export predict.rboost
#'
#' @examples
predict.rboost<- function(object,
                           newx=NULL,
                           ...) {
  if (object$rc){
    object$tau.const + predict(object$tau.fit, newx=newx)
  }
  else{
    predict(object$tau.fit, newx=newx)
  }
}
