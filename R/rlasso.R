#' R-learner, as proposed by Nie and Wager 2017, implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds number of folds for cross-fitting
#' @param lambda.choice how to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param rs whether to use the RS-learner (logical).
#' @param p_hat user-supplied estimate for E[W|X]
#' @param m_hat user-supplied estimte for E[Y|X]
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rlasso.fit = rlasso(x, w, y)
#' rlasso.est = predict(rlasso.fit, x)
#' }
#' @export
rlasso = function(X, W, Y,
                  alpha = 1,
                  k_folds=NULL,
                  lambda.choice=c("lambda.min","lambda.1se"),
                  rs = FALSE,
                  p_hat = NULL,
                  m_hat = NULL){

    if (is.null(colnames(X))) {
      stop("The design matrix X must have named columns.")
    }
    standardization = caret::preProcess(X, method=c("center", "scale")) # get the standardization params
    X.scl = predict(standardization, X)							 # standardize the input
    X.scl = X.scl[,!is.na(colSums(X.scl))]

    lambda.choice = match.arg(lambda.choice)

    nobs = nrow(X.scl)
    pobs = ncol(X.scl)

    if (is.null(k_folds)) {
        k_folds = floor(max(3, min(10,length(W)/4)))
    }

    # fold ID for cross-validation; balance treatment assignments
    foldid = sample(rep(seq(k_folds), length = length(W)))

    if (is.null(m_hat)){
      y.fit = glmnet::cv.glmnet(X, Y, foldid = foldid, keep = TRUE, alpha = alpha)
      m_hat = y.fit$fit.preval[,!is.na(colSums(y.fit$fit.preval))][, y.fit$lambda == y.fit$lambda.min]
    }
    else {
      y.fit = NULL
    }

    if (is.null(p_hat)){
      w.fit = glmnet::cv.glmnet(X, W, foldid = foldid, keep = TRUE, family = "binomial", type.measure = "deviance", alpha = alpha)
      p_hat = w.fit$fit.preval[,!is.na(colSums(w.fit$fit.preval))][, w.fit$lambda == w.fit$lambda.min]
    }
    else{
      w.fit = NULL
    }

    Y.tilde = Y - m_hat

    if (rs){

      X.scl.tilde = cbind(as.numeric(W - p_hat) * cbind(1, X.scl), X.scl)
      X.scl.pred = cbind(1, X.scl, X.scl * 0)
      penalty.factor = c(0, rep(1, 2 * pobs))

    }
    else{

      X.scl.tilde = cbind(as.numeric(W - p_hat) * cbind(1, X.scl))
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
               p_hat = p_hat,
               m_hat = m_hat,
               tau.hat = tau.hat,
               rs = rs,
               standardization = standardization)
    class(ret) <- "rlasso"
    ret
}


#' predict for rlasso
#'
#' get estimated tau(x) using the trained rlasso model
#'
#' @param object a rlasso object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param ... additional arguments (currently not used)
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rlasso.fit = rlasso(x, w, y)
#' rlasso.est = predict(rlasso.fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
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
