#' S-learner, as proposed by Imai and Ratkovic 2013
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the treatment variable (0 or 1)
#' @param nfolds number of folds for cross-fitting
#'
#' @export sboost
sboost = function(X, Y, W,
                  nfolds = NULL){

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,nobs/4)))
  }

  s.fit = cvboost(cbind(X, (W-0.5)*X, (W-0.5)), Y, objective="reg:linear", nfolds=nfolds)

  mu0.hat = predict(s.fit, newx=cbind(X, (0-0.5)*X, (0-0.5)))
  mu1.hat = predict(s.fit, newx=cbind(X, (1-0.5)*X, (1-0.5)))
  tau.hat = mu1.hat - mu0.hat

  ret = list(s.fit = s.fit,
             mu0.hat = mu0.hat,
             mu1.hat = mu1.hat,
             tau.hat = tau.hat)

  class(ret) <- "sboost"
  ret
}

#' @export predict.sboost
predict.sboost <- function(object,
                           newx=NULL,
                           ...) {
  if (!is.null(newx)) {
    mu0.hat = predict(object$s.fit, newx=cbind(newx, (0-0.5)*newx, (0-0.5)))
    mu1.hat = predict(object$s.fit, newx=cbind(newx, (1-0.5)*newx, (1-0.5)))
    tau.hat = mu1.hat - mu0.hat
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
