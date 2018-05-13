#' R-lasso, as proposed by Nie and Wager 2017
#'
#' @param X
#' @param Y
#' @param W
#' @param alpha
#' @param nfolds
#' @param lambda.choice
#' @param constant.effect
#' @param rs
#'
#' @return
#' @export rlasso
#'
#' @examples
rlasso = function(X, Y, W,
                  alpha = 1,
                  nfolds=NULL,
                  lambda.choice=c("lambda.min","lambda.1se"),
                  constant.effect = TRUE,
                  rs = FALSE,
                  w.hat = NULL,
                  y.hat = NULL) {

    X.scl = scale(X)
    X.scl = X.scl[,!is.na(colSums(X.scl))]

    lambda.choice = match.arg(lambda.choice)

    nobs = nrow(X.scl)
    pobs = ncol(X.scl)

    if (is.null(nfolds)) {
        nfolds = floor(max(3, min(10,length(W)/4)))
    }

    # fold ID for cross-validation; balance treatment assignments
    foldid = sample(rep(seq(nfolds), length = length(W)))

    if (is.null(y.hat)){
      y.fit = glmnet::cv.glmnet(X, Y, foldid=foldid, keep=TRUE, alpha = alpha)
      y.hat = y.fit$fit.preval[,!is.na(colSums(y.fit$fit.preval))][, y.fit$lambda == y.fit$lambda.min]
    }

    if (is.null(w.hat)){
      w.fit = glmnet::cv.glmnet(X, W, foldid=foldid, keep=TRUE, family="binomial", type.measure = "auc", alpha = alpha)
      w.hat = w.fit$fit.preval[,!is.na(colSums(w.fit$fit.preval))][, w.fit$lambda == w.fit$lambda.min]
    }

    Y.tilde = Y - y.hat

    if (rs){

      if (constant.effect){
        X.scl.tilde = cbind(as.numeric(W - w.hat) * cbind(1, X.scl), X.scl)
        X.scl.pred = cbind(1, X.scl, X.scl * 0)
        penalty.factor = c(0, rep(1, 2 * pobs))
      }
      else{
        X.scl.tilde = cbind(as.numeric(W - w.hat) * X.scl, X.scl)
        X.scl.pred = cbind(X.scl, X.scl * 0)
        penalty.factor = rep(1, 2 * pobs)
      }

    }
    else{

      if (constant.effect){
        X.scl.tilde = cbind(as.numeric(W - w.hat) * cbind(1, X.scl))
        X.scl.pred = cbind(1, X.scl)
        penalty.factor = c(0, rep(1, pobs))
      }
      else{
        X.scl.tilde = cbind(as.numeric(W - w.hat) * X.scl)
        X.scl.pred =  X.scl
        penalty.factor = rep(1, pobs)
      }
    }



    tau.fit = glmnet::cv.glmnet(X.scl.tilde, Y.tilde, foldid = foldid,
                             alpha = alpha,
                             penalty.factor = penalty.factor,
                             standardize = FALSE)

    tau.beta = as.vector(t(coef(tau.fit, s=lambda.choice)[-1]))

    tau.hat = X.scl.pred %*% tau.beta

    return(list(tau.hat = tau.hat, y.hat = y.hat, w.hat = w.hat, tau.beta = tau.beta))
}
