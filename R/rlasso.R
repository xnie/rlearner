#' R-lasso, as proposed by Nie and Wager 2017
#'
#' @param X
#' @param Y
#' @param W
#' @param alpha
#' @param nfolds
#' @param lambda.choice
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
                  rs = FALSE,
                  w.hat = NULL,
                  y.hat = NULL,
                  penalty.search=FALSE,
                  w.measure=c("deviance","auc"),
                  pilot.lambda.choice=c("lambda.min","lambda.1se")){

    X.scl = scale(X)
    X.scl = X.scl[,!is.na(colSums(X.scl))]

    lambda.choice = match.arg(lambda.choice)
    w.measure = match.arg(w.measure)
    pilot.lambda.measure = match.arg(pilot.lambda.choice)

    nobs = nrow(X.scl)
    pobs = ncol(X.scl)

    if (is.null(nfolds)) {
        nfolds = floor(max(3, min(10,length(W)/4)))
    }

    # fold ID for cross-validation; balance treatment assignments
    foldid = sample(rep(seq(nfolds), length = length(W)))

    if (is.null(y.hat)){
      y.fit = glmnet::cv.glmnet(X, Y, foldid=foldid, keep=TRUE, alpha = alpha)
      if (pilot.lambda.choice == "lambda.min"){
        y.hat = y.fit$fit.preval[,!is.na(colSums(y.fit$fit.preval))][, y.fit$lambda == y.fit$lambda.min]
      }
      else{
        y.hat = y.fit$fit.preval[,!is.na(colSums(y.fit$fit.preval))][, y.fit$lambda == y.fit$lambda.1se]
      }
    }
    else {
      y.fit = NULL
    }

    if (is.null(w.hat)){
      w.fit = glmnet::cv.glmnet(X, W, foldid=foldid, keep=TRUE, family="binomial", type.measure = w.measure, alpha = alpha)
      if (pilot.lambda.choice == "lambda.min"){
        w.hat = w.fit$fit.preval[,!is.na(colSums(w.fit$fit.preval))][, w.fit$lambda == w.fit$lambda.min]
      }
      else{
        w.hat = w.fit$fit.preval[,!is.na(colSums(w.fit$fit.preval))][, w.fit$lambda == w.fit$lambda.1se]
      }
    }
    else{
      w.fit = NULL
    }

    Y.tilde = Y - y.hat

    if (rs){

      X.scl.tilde = cbind(as.numeric(W - w.hat) * cbind(1, X.scl), X.scl)
      X.scl.pred = cbind(1, X.scl, X.scl * 0)

      if(penalty.search){
        search.range = 5
        cvm.min = Inf
        last.best = NULL
        for (l in 0:2){
          updated = FALSE

          for (i in 1:search.range){

            if (l==0){
              penalty.factor = c(0, rep(10^(i- ceiling(search.range/2.0)), pobs), rep(1,pobs))
            }
            else{
              penalty.factor = c(0, rep(10^(last.best - 10^(-l+1) + 20.0/search.range*10^(-l)*i), pobs), rep(1,pobs))
            }
            rs.fit <- glmnet::cv.glmnet(x=X.scl.tilde, y=Y, foldid=foldid, penalty.factor=penalty.factor, standardize=FALSE)
            rs.fit.cvm = rs.fit$cvm[rs.fit$lambda == rs.fit$lambda.min]
            if (rs.fit.cvm < cvm.min){
              cvm.min = rs.fit.cvm
              rs.fit.best = rs.fit
              if (l==0){
                best.i = i - ceiling(search.range/2.0)
              }
              else{
                best.i = i
              }
              best.penalty.factor= penalty.factor
              updated = TRUE
            }
          }

          if (l==0){
            last.best = best.i
          }
          else{
            if (updated){
              last.best = last.best - 10^(-l+1) + 20.0/search.range*10^(-l)*best.i
            }
          }
        }
        penalty.factor = best.penalty.factor
      }
      else{
        penalty.factor = c(0, rep(1, 2 * pobs))
      }
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
               rs = rs)
    class(ret) <- "rlasso"
    ret

}

predict.rlasso <- function(object,
                           newx=NULL,
                           ...) {
  newx.scl = scale(newx)
  newx.scl = newx.scl[,!is.na(colSums(newx.scl))]
  if (!is.null(newx)) {
    if (object$rs){
      newx.scl.pred = cbind(1, newx.scl, newx.scl * 0)
    }
    else{
      newx.scl.pred = cbind(1, newx.scl)
    }
    tau.hat = newx.scl.pred %*% object$tau.beta # TODO: make this more robust?
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
