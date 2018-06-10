#' Title
#'
#' @param X
#' @param Y
#' @param weights
#' @param nfolds
#' @param objective
#' @param ntrees.max
#' @param num.search.rounds
#' @param print.every.n
#' @param early.stopping.rounds
#' @param bayes.opt
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
                   ntrees.max=1000,
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


  dtrain <- xgboost::xgb.DMatrix(data = X, label = Y, weight = weights)

  best.param = list()
  best.seednumber = 1234
  best.loss = Inf


  if (bayes.opt){ # WARNING: very slow!
    xgb_cv_bayes <- function(subsample, eta, max_depth, min_child_weight) {
      cv <- xgboost::xgb.cv(params = list(subsample = subsample,
                                          eta = eta,
                                          max_depth = max_depth,
                                          min_child_weight = min_child_weight,
                                          lambda = 1,
                                          alpha = 0,
                                          objective = objective,
                                          eval_metric = eval),
                            data = dtrain,
                            nround = ntrees.max,
                            #nthread = 1,
                            nfold = nfolds,
                            prediction = TRUE,
                            showsd = TRUE,
                            early_stopping_rounds = 10,
                            maximize = FALSE,
                            print_every_n = print.every.n,
                            callbacks = list(cb.cv.predict(save_models = TRUE)))

      metric = paste('test_', eval, '_mean', sep='')
      list(Score = min(cv$evaluation_log[, ..metric]),
           Pred = cv$pred)
    }
    opt.res <- rBayesianOptimization::BayesianOptimization(xgb_cv_bayes,
                                                           bounds = list(subsample = c(0.5, 1),
                                                                         eta = c(5e-3, 2e-1),
                                                                         max_depth = c(3L, 20L),
                                                                         min_child_weight = c(1L, 20L)),
                                                           init_grid_dt = NULL, init_points = 10, n_iter = 1,
                                                           acq = "ucb", kappa = 2.576, eps = 0.0,
                                                           verbose = TRUE)
    best.param = as.list(opt.res$Best_Par)
    best.loss = opt.res$Best_Value
    best.xgb.cvfit <- xgboost::xgb.cv(params = best.param,
                                      data = dtrain,
                                      nround = ntrees.max,
                                      #nthread = 1,
                                      nfold = nfolds,
                                      prediction = TRUE,
                                      showsd = TRUE,
                                      early_stopping_rounds = 10,
                                      maximize = FALSE,
                                      print_every_n = print.every.n,
                                      callbacks = list(cb.cv.predict(save_models = TRUE)))

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
      xgb.cvfit <- xgboost::xgb.cv(data = dtrain,
                                   param =param,
                                   missing = NA,
                                   nfold = nfolds,
                                   nthread=1,
                                   prediction = TRUE,
                                   early_stopping_rounds = early.stopping.rounds,
                                   maximize = FALSE,
                                   nrounds = ntrees.max,
                                   print_every_n = print.every.n,
                                   callbacks = list(cb.cv.predict(save_models = TRUE)))

      metric = paste('test_', eval, '_mean', sep='')

      min.loss = min(xgb.cvfit$evaluation_log[, metric])

      if (min.loss < best.loss) {
        best.loss = min.loss
        best.seednumber = seed.number
        best.param = param
        best.xgb.cvfit = xgb.cvfit
      }
    }
  }


 set.seed(best.seednumber)
 xgb.fit <- xgboost::xgb.train(data=dtrain, params=best.param, nrounds=best.xgb.cvfit$best_ntreelimit, nthread=1)

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
   dtest <- xgboost::xgb.DMatrix(data=newx)
   return(predict(object$xgb.fit, newdata=dtest))
  }
}
