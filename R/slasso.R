#' S-lasso
#'
#' @param X
#' @param Y
#' @param W
#' @param alpha
#' @param nfolds
#' @param lambda.choice
#' @param standardize
#' @param constant.effect
#'
#' @return
#' @export slasso
#'
#' @examples
slasso = function(X, Y, W,
                  alpha = 1,
                  nfolds=NULL,
                  lambda.choice=c("lambda.1se", "lambda.min"),
                  standardize = FALSE,
                  constant.effect = TRUE) {
  # standardize in glmnet
  if (standardize){
    X.scl = X
  }
  else{
    X.scl = scale(X)
    X.scl = X.scl[,!is.na(colSums(X.scl))]
  }

  lambda.choice = match.arg(lambda.choice)

  nobs = nrow(X.scl)
  pobs = ncol(X.scl)

  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,nobs/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid = sample(rep(seq(nfolds), length = nobs))

  if (constant.effect){
    X.scl.tilde = cbind(as.numeric(2 * W - 1) * cbind(1, X.scl), X.scl)
    X.scl.pred = cbind(1, X.scl, 0 * X.scl)
    penalty.factor = c(0, rep(1, 2 * pobs))
  }
  else{
    X.scl.tilde = cbind(as.numeric(2 * W - 1) * X.scl, X.scl)
    X.scl.pred = cbind(X.scl, 0 * X.scl)
    penalty.factor = rep(1, 2 * pobs)
  }


  s.fit = cv.glmnet(X.scl.tilde, Y, foldid = foldid,
                      penalty.factor = penalty.factor,
                      standardize=standardize)

  tau.hat = 2 * X.scl.pred %*% as.vector(t(coef(s.fit, s=lambda.choice)[-1]))

  return(list(tau.hat = tau.hat))
}
