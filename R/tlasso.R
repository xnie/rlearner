#' T-lasso
#'
#' @param X
#' @param Y
#' @param W
#' @param alpha
#' @param nfolds.1
#' @param nfolds.0
#' @param lambda.choice
#'
#' @return
#' @export
#'
#' @examples
tlasso = function(X, Y, W,
                  alpha = 1,
                  nfolds.1=NULL,
                  nfolds.0=NULL,
                  lambda.choice=c("lambda.min", "lambda.1se")){

  lambda.choice = match.arg(lambda.choice)

  X.1 = X[which(W==1),]
  X.0 = X[which(W==0),]

  Y.1 = Y[which(W==1)]
  Y.0 = Y[which(W==0)]

  nobs.1 = nrow(X.1)
  nobs.0 = nrow(X.0)

  pobs = ncol(X)

  if (is.null(nfolds.1)) {
    nfolds.1 = floor(max(3, min(10,nobs.1/4)))
  }

  if (is.null(nfolds.0)) {
    nfolds.0 = floor(max(3, min(10,nobs.0/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid.1 = sample(rep(seq(nfolds.1), length = nobs.1))
  foldid.0 = sample(rep(seq(nfolds.0), length = nobs.0))

  t.1.fit = glmnet::cv.glmnet(X.1, Y.1, foldid = foldid.1)
  t.0.fit = glmnet::cv.glmnet(X.0, Y.0, foldid = foldid.0)

  y.1.pred = predict(t.1.fit, newx=X, s=lambda.choice)
  y.0.pred = predict(t.0.fit, newx=X, s=lambda.choice)

  tau.hat = y.1.pred - y.0.pred

  return(list(tau.hat = tau.hat))
}
