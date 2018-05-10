#' Estimate treatment effect via the R-lasso, as proposed by Nie and Wager (2017)
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the effect variable (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param nfolds number of folds for cross-fitting
#' @param lambda.choice how to cross-validated
#' @param standardize whether X should be rescaled before running the lasso
#'
#' @export rlasso
rlasso = function(X, Y, W,
                  alpha = 1,
                  nfolds=NULL,
                  lambda.choice=c("lambda.1se", "lambda.min"),
                  standardize = TRUE) {
    
    lambda.choice = match.arg(lambda.choice)
    
    nobs = nrow(X)
    pobs = ncol(X)
    
    if (is.null(nfolds)) {
        nfolds = floor(max(3, min(10,length(W)/4)))
    }
    
    # fold ID for cross-validation; balance treatment assignments
    foldid = sample(rep(seq(nfolds), length = length(W)))
    
    y.fit = cv.glmnet(X, Y, foldid=foldid, keep=TRUE, alpha = alpha, standardize = standardize)
    y.hat = y.fit$fit.preval[, y.fit$lambda == y.fit$lambda.min]
    
    w.fit = cv.glmnet(X, W, foldid=foldid, keep=TRUE, family="binomial", type.measure = "auc", alpha = alpha, standardize = standardize)
    w.hat = w.fit$fit.preval[, w.fit$lambda == w.fit$lambda.min]
    
    Y.tilde = Y - y.hat
    X.tilde = cbind(as.numeric(W - w.hat) * cbind(1, X)) # TODO: how to add lambda choice to cv glmnet?

    tau.fit = cv.glmnet(X.tilde, Y.tilde, foldid = foldid,
                             alpha = alpha,
                             penalty.factor = c(0, rep(1, pobs)),
                             standardize = standardize)
    
    tau.beta = as.vector(t(coef(tau.fit, s=lambda.choice)[-1]))
    tau.hat = cbind(1, X) %*% tau.beta
    print(tau.hat)
    
    return(list(tau.hat = tau.hat, y.hat = y.hat, w.hat = w.hat, tau.beta = tau.beta))
}
