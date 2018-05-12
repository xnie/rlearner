#' U-lasso, as proposed by Kunzel et al 2017
#'
#' @param X
#' @param Y
#' @param W
#' @param alpha
#' @param nfolds
#' @param lambda.choice
#' @param cutoff
#'
#' @return
#' @export ulasso
#'
#' @examples
ulasso = function(X, Y, W,
                  alpha = 1,
                  nfolds=NULL,
                  lambda.choice=c("lambda.min", "lambda.1se"),
                  cutoff=0){

  lambda.choice = match.arg(lambda.choice)

  nobs = nrow(X)
  pobs = ncol(X)

  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,length(W)/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid = sample(rep(seq(nfolds), length = length(W)))

  y.fit = glmnet::cv.glmnet(X, Y, foldid=foldid, keep=TRUE, alpha = alpha)
  y.hat = y.fit$fit.preval[, y.fit$lambda == y.fit$lambda.min]

  w.fit = glmnet::cv.glmnet(X, W, foldid=foldid, keep=TRUE, family="binomial", type.measure = "auc", alpha = alpha)
  w.hat = w.fit$fit.preval[, w.fit$lambda == w.fit$lambda.min]
  w.hat.thresh = pmax(cutoff, pmin(1 - cutoff, w.hat))

  Y.tilde = Y - y.hat
  W.tilde = W - w.hat.thresh

  U = Y.tilde / W.tilde

  tau.fit = glmnet::cv.glmnet(X, U, foldid = foldid, alpha = alpha)
  tau.hat = predict(tau.fit, newx = X)

  return(list(tau.hat = tau.hat))
}
