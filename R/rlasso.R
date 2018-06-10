#' R-learner, as proposed by Nie and Wager, 2017
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the treatment variable (0 or 1)
#' @param alpha tuning parameter for the elastic net
#' @param nfolds number of folds for cross-fitting
#' @param lambda.choice how to cross-validate
#' @param rs whether to use the RS-learner (logical)
#' @param w.hat user-supplied estimate for E[W|X]
#' @param y.hat user-supplied estimte for E[Y|X]
#'
#' @export rlasso
rlasso = function(X, Y, W,
                  alpha = 1,
                  nfolds=NULL,
                  lambda.choice=c("lambda.min","lambda.1se"),
                  rs = FALSE,
                  w.hat = NULL,
                  y.hat = NULL){

    #X.scl = scale(X)
    #X.scl = X.scl[,!is.na(colSums(X.scl))]

    if (is.null(colnames(X))) {
      stop("The design matrix X must have named columns.")
    }
    standardization = caret::preProcess(X, method=c("center", "scale")) # get the standardization params
    X.scl = predict(standardization, X)							 # standardize the input
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
      y.fit = glmnet::cv.glmnet(X, Y, foldid = foldid, keep = TRUE, alpha = alpha)
      y.hat = y.fit$fit.preval[,!is.na(colSums(y.fit$fit.preval))][, y.fit$lambda == y.fit$lambda.min]
    }
    else {
      y.fit = NULL
    }

    if (is.null(w.hat)){
      w.fit = glmnet::cv.glmnet(X, W, foldid = foldid, keep = TRUE, family = "binomial", type.measure = "deviance", alpha = alpha)
      w.hat = w.fit$fit.preval[,!is.na(colSums(w.fit$fit.preval))][, w.fit$lambda == w.fit$lambda.min]
    }
    else{
      w.fit = NULL
    }

    Y.tilde = Y - y.hat

    if (rs){

      X.scl.tilde = cbind(as.numeric(W - w.hat) * cbind(1, X.scl), X.scl)
      X.scl.pred = cbind(1, X.scl, X.scl * 0)
      penalty.factor = c(0, rep(1, 2 * pobs))

    }
    else{

      X.scl.tilde = cbind(as.numeric(W - w.hat) * cbind(1, X.scl))
      X.scl.pred = cbind(1, X.scl)
      penalty.factor = c(0, rep(1, pobs))

    }

    tau.fit = glmnet::cv.glmnet(X.scl.tilde, Y.tilde, foldid = foldid,
                             alpha = alpha,
                             penalty.factor = penalty.factor,
                             standardize = FALSE)

    tau.beta = as.vector(t(coef(tau.fit, s=lambda.choice)[-1]))

    tau.hat = X.scl.pred %*% tau.beta

    ret = list(tau.fit = tau.fit,
               tau.beta = tau.beta,
               w.fit = w.fit,
               y.fit = y.fit,
               w.hat = w.hat,
               y.hat = y.hat,
               tau.hat = tau.hat,
               rs = rs,
               standardization = standardization)
    class(ret) <- "rlasso"
    ret
}

#' Title
#'
#' @param object
#' @param newx
#' @param ...
#'
#' @return
#' @export predict.rlasso
#'
#' @examples
predict.rlasso <- function(object,
                           newx=NULL,
                           ...) {
  if (!is.null(newx)) {
    newx.scl = predict(object$standardization, newx) # standardize the new data using the same standardization as with the training data
    newx.scl = newx.scl[,!is.na(colSums(newx.scl))]

    if (object$rs){
      newx.scl.pred = cbind(1, newx.scl, newx.scl * 0)
    }
    else{
      newx.scl.pred = cbind(1, newx.scl)
    }
    tau.hat = newx.scl.pred %*% object$tau.beta
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
