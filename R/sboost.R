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

  nobs = nrow(X.scl)
  pobs = ncol(X.scl)

  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,nobs/4)))
  }

  browser()

  s.fit = cvboost(cbind(X, (W-0.5)*X, (W-0.5)), Y, objective="reg:linear")

  list(0, 1) %>% purrr::map(function(w) {
    predict(s.fit, newx=cbind(X, (w-0.5)*X, (W-0.5)))
  }) %->% c(mu0.hat, mu1.hat)
  tau.hat = mu1.hat - mu0.hat

  ret = list(s.fit = s.fit,
             tau.hat = tau.hat)

  class(ret) <- "sboost"
  ret
}

#' @export predict.sboost
predict.sboost <- function(object,
                           newx=NULL,
                           ...) {
  if (!is.null(newx)) {
    list(0, 1) %>% purrr::map(function(w) {
      predict(object$s.fit, newdata=cbind(newx, (w-0.5)*newx, (w-0.5)))
    }) %->% c(mu0.hat, mu1.hat)
    tau.hat = mu1.hat - mu0.hat
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
