#require(xgboost)
#X = matrix(rnorm(5000), 500, 10)
#b = pmax(0, X[,1] + X[,2], X[,3]) + pmax(0, X[,4] + X[,5])
#e = 0.5
#tau = X[,1] + log(1 + exp(X[,2]))
#W = rbinom(500,1,e)==1
#Y = b + (W - 0.5) * tau + rnorm(500)
#weights=NULL
#nfolds=NULL
#objective="reg:linear"
#ntrees.max=1000

#' Title
#'
#' @param X
#' @param Y
#' @param weights
#' @param nfolds
#' @param objective
#' @param ntrees.max
#' @param num.search.rounds
#'
#' @return
#' @export
#'
#' @examples
cvboost = function(X,
                   Y,
                   weights=NULL,
                   nfolds=NULL,
                   objective=c("reg:linear", "binary:logistic"),
                   ntrees.max=30000,
                   num.search.rounds=10,
                   print.every.n=100,
                   early.stopping.rounds=10,
                   bayes.opt=FALSE) {

  objective = match.arg(objective)
  if (objective == "reg:linear") {
    eval = "rmse"
  }
  else if (objective == "binary:logistic") {
    eval = "logloss"
  }
  else {
    stop("objective not defined.")
  }

  if (is.null(nfolds)) {
    nfolds = floor(max(3, min(10,length(Y)/4)))
  }
  if (is.null(weights)) {
    weights = rep(1, length(Y))
  }


  dtrain <- xgb.DMatrix(data = X, label = Y, weight = weights)

  best.param = list()
  best.seednumber = 1234
  best.loss = Inf
  best.loss_index = 0

  if (bayes.opt){
      #res0 <- xgb_cv_opt(data = dtrain,
      #                   label = Y,
      #                   objectfun = objective,
      #                   evalmetric = eval,
      #                   n_folds = nfolds,
      #                   eta_range = c(0.0001, 0.05),
      #                   max_depth_range = c(2L, 10L),
      #                   nrounds_range = c(10, 1000L),
      #                   subsample_range = c(0.1, 1L),
      #                   bytree_range = c(0.4, 1L),
      #                   n_iter = 10)
    stop("Bayesian optimization still needs to be added.")

  }
  else{
    for (iter in 1:num.search.rounds) {
      param <- list(objective = objective,
                    eval_metric = eval,
                    subsample = sample(c(0.5, 0.75, 1), 1),
                    colsample_bytree = sample(c(0.6, 0.8, 1), 1),
                    eta = sample(c(5e-3, 1e-2, 0.015, 0.025, 5e-2, 8e-2, 1e-1, 2e-1), 1),
                    max_depth = sample(c(3:20), 1),
                    gamma = runif(1, 0.0, 0.2),
                    min_child_weight = sample(1:20, 1),
                    max_delta_step = sample(1:10, 1))

      seed.number = sample.int(100000, 1)[[1]]
      set.seed(seed.number)
      xgb.cvfit <- xgb.cv(data = dtrain,
                          param =param,
                          missing = NA,
                          nfold = nfolds,
                          prediction = TRUE,
                          early_stopping_rounds = early.stopping.rounds,
                          maximize = FALSE,
                          nrounds = ntrees.max,
                          print_every_n = print.every.n,
                          callbacks = list(cb.cv.predict(save_models = TRUE)))

      metric = paste('test_', eval, '_mean', sep='')

      min.loss = min(xgb.cvfit$evaluation_log[, metric])
      #min.loss.index = xgb.cvfit$best_iteration

      if (min.loss < best.loss) {
        best.loss = min.loss
        #best.loss.index = min.loss.index
        best.seednumber = seed.number
        best.param = param
        best.xgb.cvfit = xgb.cvfit
      }
    }
  }


 set.seed(best.seednumber)
 xgb.fit <- xgb.train(data=dtrain, params=best.param, nrounds=best.xgb.cvfit$best_ntreelimit)

 ret = list(xgb.fit = xgb.fit,
            best.xgb.cvfit = best.xgb.cvfit,
            best.seednumber = best.seednumber,
            best.param = best.param,
            best.loss = best.loss,
            best.ntreelimit = best.xgb.cvfit$best_ntreelimit)

  class(ret) <- "cvboost"
  ret
}

#' Title
#'
#' @param object
#' @param newx
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict.cvboost <- function(object,
                            newx=NULL,
                            ...) {
  if (is.null(newx)) {
    return(object$best.xgb.cvfit$pred)
  }
  else{
   dtest <- xgb.DMatrix(data=newx)
   return(predict(object$xgb.fit, newdata=dtest))
  }
}
