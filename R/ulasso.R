#' U-learner
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the treatment variable (0 or 1)
#' @param alpha tuning parameter for the elastic net
#' @param nfolds number of folds for cross-fitting
#' @param lambda.choice how to cross-validate
#' @param cutoff propensity cutoff threshold
#'
#' @export ulasso
ulasso = function(X, Y, W,
                  alpha = 1,
                  nfolds=NULL,
                  lambda.choice=c("lambda.1se", "lambda.mse"),
                  cutoff=0.05){

  lambda.choice = match.arg(lambda.choice)

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,length(W)/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid = sample(rep(seq(nfolds), length = length(W)))

  y.fit = glmnet::cv.glmnet(X, Y, nfolds=10, keep=TRUE, alpha = alpha)
  y.hat = y.fit$fit.preval[,!is.na(colSums(y.fit$fit.preval))][, y.fit$lambda == y.fit$lambda.min]

  w.fit = glmnet::cv.glmnet(X, W, nfolds=10, keep=TRUE, family="binomial", type.measure = "deviance", alpha = alpha)
  w.hat = w.fit$fit.preval[,!is.na(colSums(w.fit$fit.preval))][, w.fit$lambda == w.fit$lambda.min]

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

#' @export predict.ulasso
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
