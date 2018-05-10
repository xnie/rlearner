#' R-lasso, as proposed by Nie and Wager 2017
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
#' @export rlasso
#'
#' @examples
rlasso = function(X, Y, W,
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
        nfolds = floor(max(3, min(10,length(W)/4)))
    }

    # fold ID for cross-validation; balance treatment assignments
    foldid = sample(rep(seq(nfolds), length = length(W)))

    y.fit = cv.glmnet(X, Y, foldid=foldid, keep=TRUE, alpha = alpha)
    y.hat = y.fit$fit.preval[, y.fit$lambda == y.fit$lambda.min]

    w.fit = cv.glmnet(X, W, foldid=foldid, keep=TRUE, family="binomial", type.measure = "auc", alpha = alpha)
    w.hat = w.fit$fit.preval[, w.fit$lambda == w.fit$lambda.min]

    Y.tilde = Y - y.hat

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

    tau.fit = cv.glmnet(X.scl.tilde, Y.tilde, foldid = foldid,
                             alpha = alpha,
                             penalty.factor = penalty.factor,
                             standardize = standardize)

    tau.beta = as.vector(t(coef(tau.fit, s=lambda.choice)[-1]))
    tau.hat = X.scl.pred %*% tau.beta

    return(list(tau.hat = tau.hat, y.hat = y.hat, w.hat = w.hat, tau.beta = tau.beta))
}
