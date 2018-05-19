#' S-learner, as proposed by Imai and Ratkovic 2013
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the treatment variable (0 or 1)
#' @param alpha tuning parameter for the elastic net
#' @param nfolds number of folds for cross-fitting
#' @param lambda.choice how to cross-validate
#' @param penalty.search whether to perform fine grainted penalty factor search (logical)
#'
#' @export slasso
slasso = function(X, Y, W,
                  alpha = 1,
                  nfolds = NULL,
                  lambda.choice = c("lambda.min", "lambda.1se"),
                  penalty.search = FALSE) {

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
  if(penalty.search){
    search.range = 5
    cvm.min = Inf
    last.best = NULL
    for (l in 0:2){
      #print(paste("level :",l))

      updated = FALSE

      for (i in 1:search.range){

        if (l==0){
          penalty.factor = c(0, rep(10^(i- ceiling(search.range/2.0)), pobs), rep(1,pobs))
        }
        else{
          penalty.factor = c(0, rep(10^(last.best - 10^(-l+1) + 20.0/search.range*10^(-l)*i), pobs), rep(1,pobs))
        }
        s.fit <- glmnet::cv.glmnet(x = X.scl.tilde, y = Y, foldid = foldid, penalty.factor = penalty.factor, standardize = FALSE, alpha = alpha)
        s.fit.cvm = s.fit$cvm[s.fit$lambda == s.fit$lambda.min]
        s.beta = as.vector(t(coef(s.fit, s=lambda.choice)[-1]))
        if (s.fit.cvm < cvm.min){
          cvm.min = s.fit.cvm
          s.fit.best = s.fit
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

  s.fit = glmnet::cv.glmnet(X.scl.tilde, Y, foldid = foldid,
                      penalty.factor = penalty.factor,
                      standardize = FALSE, alpha = alpha)

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
