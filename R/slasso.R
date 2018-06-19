#' S-learner, as proposed by Imai and Ratkovic 2013, implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds number of folds for cross validation
#' @param lambda.choice how to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty.search whether to perform fine grainted penalty factor search (logical)
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' slasso.fit = slasso(x, w, y)
#' slasso.est = predict(slasso.fit, x)
#' }
#' @export
slasso = function(X, W, Y,
                  alpha = 1,
                  k_folds = NULL,
                  lambda.choice = c("lambda.min", "lambda.1se"),
                  penalty.search = FALSE) {

  standardization = caret::preProcess(X, method=c("center", "scale")) # get the standardization params
  X.scl = predict(standardization, X)							 # standardize the input
  X.scl = X.scl[,!is.na(colSums(X.scl))]

  lambda.choice = match.arg(lambda.choice)

  nobs = nrow(X.scl)
  pobs = ncol(X.scl)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,nobs/4)))
  }

  # fold ID for cross-validation; balance treatment assignments
  foldid = sample(rep(seq(k_folds), length = nobs))

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
             tau.hat = tau.hat,
             standardization = standardization)


  class(ret) <- "slasso"
  ret
}

#' predict for slasso
#'
#' get estimated tau(x) using the trained slasso model
#'
#' @param object a slasso object
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
#' slasso.fit = slasso(x, w, y)
#' slasso.est = predict(slasso.fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.slasso <- function(object,
                           newx=NULL,
                           ...) {
  if (!is.null(newx)) {
    newx.scl = predict(object$standardization, newx) # standardize the new data using the same standardization as with the training data
    newx.scl = newx.scl[,!is.na(colSums(newx.scl))]
    newx.scl.pred = cbind(1, newx.scl, 0 * newx.scl)
    tau.hat = 2 * newx.scl.pred %*% object$s.beta
  }
  else {
    tau.hat = object$tau.hat
  }
  return(tau.hat)
}
