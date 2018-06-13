#' Title
#'
#' @param X
#' @param Y
#' @param W
#' @param nfolds
#' @param w.hat
#' @param y.hat
#' @param cutoff
#'
#' @return
#' @export
#'
#' @examples
uboost= function(X, Y, W,
                 nfolds=NULL,
                 w.hat = NULL,
                 y.hat = NULL,
                 cutoff=0.05,
                 nthread=NULL){

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,length(W)/4)))
  }

  if (is.null(y.hat)){
    y.fit = cvboost(X, Y, objective="reg:linear", nfolds=nfolds, nthread=nthread)
    y.hat = predict(y.fit)
  }
  else {
    y.fit = NULL
  }

  if (is.null(w.hat)){
    w.fit = cvboost(X, W, objective="binary:logistic", nfolds=nfolds, nthread=nthread)
    w.hat = predict(w.fit)
  }
  else{
    w.fit = NULL
  }

  w.hat = pmax(cutoff, pmin(1 - cutoff, w.hat))

  Y.tilde = Y - y.hat
  W.tilde = W - w.hat
  pseudo.outcome = Y.tilde/W.tilde
  tau.const = NULL


  tau.fit = cvboost(X, pseudo.outcome, objective="reg:linear", nfolds=nfolds, nthread=nthread)

  ret = list(tau.fit = tau.fit,
             w.fit = w.fit,
             y.fit = y.fit,
             w.hat = w.hat,
             y.hat = y.hat,
             tau.const = tau.const)
  class(ret) <- "uboost"
  ret
}

#' Title
#'
#' @param object
#' @param newx
#' @param ...
#'
#' @return
#' @export predict.uboost
#'
#' @examples
predict.uboost<- function(object,
                          newx=NULL,
                          ...) {
  predict(object$tau.fit, newx=newx)
}
