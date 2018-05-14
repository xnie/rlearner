#' S-lasso
#'
#' @param X
#' @param Y
#' @param W
#' @param alpha
#' @param nfolds
#' @param lambda.choice
#' @param constant.effect
#'
#' @return
#' @export slasso
#'
#' @examples
slasso = function(X, Y, W,
                  alpha = 1,
                  nfolds=NULL,
                  lambda.choice=c("lambda.min", "lambda.1se")) {

  X.scl = scale(X)
  X.scl = X.scl[,!is.na(colSums(X.scl))]

  lambda.choice = match.arg(lambda.choice)

  nobs = nrow(X.scl)
  pobs = ncol(X.scl)

  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,nobs/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid = sample(rep(seq(nfolds), length = nobs))

  X.scl.tilde = cbind(as.numeric(2 * W - 1) * cbind(1, X.scl), X.scl)
  X.scl.pred = cbind(1, X.scl, 0 * X.scl)
  penalty.factor = c(0, rep(1, 2 * pobs))

  s.fit = glmnet::cv.glmnet(X.scl.tilde, Y, foldid = foldid,
                      penalty.factor = penalty.factor,
                      standardize=FALSE)

  s.beta = as.vector(t(coef(s.fit, s=lambda.choice)[-1]))

  tau.hat = 2 * X.scl.pred %*% s.beta

  ret = list(s.fit = s.fit,
             s.beta = s.beta,
             tau.hat = tau.hat)

  class(ret) <- "slasso"
  ret
}

predict.slasso <- function(object,
                           newx=NULL,
                           ...) {
  newx.scl = scale(newx)
  newx.scl = newx.scl[,!is.na(colSums(newx.scl))]
  newx.scl.pred = cbind(1, newx.scl, 0 * newx.scl)

  if (!is.null(newx)) {
    tau.hat = 2 * newx.scl.pred %*% object$s.beta
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
