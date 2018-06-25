#' Gradient boosting for regression and classification with cross validation to search for hyper-parameters (implemented with xgboost)
#'
#' @param x the input features
#' @param y the observed response (real valued)
#' @param weights weights for input if doing weighted regression/classification. If set to NULL, no weights are used
#' @param k_folds number of folds used in cross validation
#' @param objective choose from either "reg:linear" for regression or "binary:logistic" for logistic regression
#' @param ntrees_max the maximum number of trees to grow for xgboost
#' @param num_search_rounds the number of random sampling of hyperparameter combinations for cross validating on xgboost trees
#' @param print_every_n the number of iterations (in each iteration, a tree is grown) by which the code prints out information
#' @param early_stopping_rounds the number of rounds the test error stops decreasing by which the cross validation in finding the optimal number of trees stops
#' @param nthread the number of threads to use. The default is NULL, which uses all available threads. Note that this does not apply to using bayesian optimization to search for hyperparameters.
#' @param verbose boolean; whether to print statistic
#' @param bayes_opt if set to TRUE, use bayesian optimization to do hyper-parameter search in xgboost (warning: very slow!). if set to FALSE, randomly draw combinations of hyperparameters to search from (as specified by num_search_rounds). Default is FALSE.
#'
#' @return a cvboost object
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' y = pmax(x[,1], 0) + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' fit = cvboost(x, y, objective="reg:linear")
#' est = predict(fit, x)
#' }
#'
#' @export
cvboost = function(x,
                   y,
                   weights=NULL,
                   k_folds=NULL,
                   objective=c("reg:linear", "binary:logistic"),
                   ntrees_max=1000,
                   num_search_rounds=10,
                   print_every_n=100,
                   early_stopping_rounds=10,
                   nthread=NULL,
                   verbose=FALSE,
                   bayes_opt=FALSE) {

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

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(y)/4)))
  }
  if (is.null(weights)) {
    weights = rep(1, length(y))
  }

  dtrain <- xgboost::xgb.DMatrix(data = x, label = y, weight = weights)

  best_param = list()
  best_seednumber = 1234
  best_loss = Inf

  if (is.null(nthread)){
    nthread = parallel::detectCores()
  }

  if (bayes_opt){ # WARNING: in beta version; very slow!
    xgb_cv_bayes <- function(subsample, eta, max_depth, min_child_weight) {
      xgb_cv_args = list(params = list(subsample = subsample,
                                       eta = eta,
                                       max_depth = max_depth,
                                       min_child_weight = min_child_weight,
                                       lambda = 1,
                                       alpha = 0,
                                       objective = objective,
                                       eval_metric = eval),
                         data = dtrain,
                         nround = ntrees_max,
                         nfold = k_folds,
                         prediction = TRUE,
                         showsd = TRUE,
                         early_stopping_rounds = 10,
                         maximize = FALSE,
                         print_every_n = print_every_n,
                         verbose = verbose,
                         nthread = nthread,
                         callbacks = list(xgboost::cb.cv.predict(save_models = TRUE)))

      cv <- do.call(xgboost::xgb.cv, xgb_cv_args)

      metric = paste('test_', eval, '_mean', sep='')
      list(Score = min(cv$evaluation_log[, metric]),
           Pred = cv$pred)
    }
    opt_res <- rBayesianOptimization::BayesianOptimization(xgb_cv_bayes,
                                                           bounds = list(subsample = c(0.5, 1),
                                                                         eta = c(5e-3, 2e-1),
                                                                         max_depth = c(3L, 20L),
                                                                         min_child_weight = c(1L, 20L)),
                                                           init_grid_dt = NULL, init_points = 10, n_iter = num_search_rounds,
                                                           acq = "ucb", kappa = 2.576, eps = 0.0,
                                                           verbose = TRUE)
    best_param = as.list(opt_res$Best_Par)
    best_loss = opt_res$Best_Value

    xgb_cvfit_args = list(params = best_param,
                          data = dtrain,
                          nround = ntrees_max,
                          nfold = k_folds,
                          prediction = TRUE,
                          showsd = TRUE,
                          early_stopping_rounds = 10,
                          maximize = FALSE,
                          print_every_n = print_every_n,
                          verbose = verbose,
                          nthread = nthread,
                          callbacks = list(xgboost::cb.cv.predict(save_models = TRUE)))

    best_xgb_cvfit <- do.call(xgboost::xgb.cv, xgb_cvfit_args)

  }
  else{
    for (iter in 1:num_search_rounds) {
      param <- list(objective = objective,
                    eval_metric = eval,
                    subsample = sample(c(0.5, 0.75, 1), 1),
                    colsample_bytree = sample(c(0.6, 0.8, 1), 1),
                    eta = sample(c(5e-3, 1e-2, 0.015, 0.025, 5e-2, 8e-2, 1e-1, 2e-1), 1),
                    max_depth = sample(c(3:20), 1),
                    gamma = runif(1, 0.0, 0.2),
                    min_child_weight = sample(1:20, 1),
                    max_delta_step = sample(1:10, 1))

      seed_number = sample.int(100000, 1)[[1]]
      set.seed(seed_number)
      xgb_cv_args = list(data = dtrain,
                         param = param,
                         missing = NA,
                         nfold = k_folds,
                         prediction = TRUE,
                         early_stopping_rounds = early_stopping_rounds,
                         maximize = FALSE,
                         nrounds = ntrees_max,
                         print_every_n = print_every_n,
                         verbose = verbose,
                         nthread = nthread,
                         callbacks = list(xgboost::cb.cv.predict(save_models = TRUE)))

      xgb_cvfit <- do.call(xgboost::xgb.cv, xgb_cv_args)

      metric = paste('test_', eval, '_mean', sep='')

      min_loss = min(xgb_cvfit$evaluation_log[, metric])

      if (min_loss < best_loss) {
        best_loss = min_loss
        best_seednumber = seed_number
        best_param = param
        best_xgb_cvfit = xgb_cvfit
      }
    }
  }

  set.seed(best_seednumber)

  xgb_train_args = list(data = dtrain,
                        params = best_param,
                        nthread = nthread,
                        nrounds = best_xgb_cvfit$best_ntreelimit)

  xgb_fit <- do.call(xgboost::xgb.train, xgb_train_args)

  ret = list(xgb_fit = xgb_fit,
             best_xgb_cvfit = best_xgb_cvfit,
             best_seednumber = best_seednumber,
             best_param = best_param,
             best_loss = best_loss,
             best_ntreelimit = best_xgb_cvfit$best_ntreelimit)

  class(ret) <- "cvboost"
  ret
}

#' predict for cvboost
#'
#' @param object a cvboost object
#' @param newx covariate matrix to make predictions on. If null, return the predictions on the training data
#' @param ... additional arguments (currently not used)
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' y = pmax(x[,1], 0) + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' fit = cvboost(x, y, objective="reg:linear")
#' est = predict(fit, x)
#' }
#'
#' @return vector of predictions
#' @export
predict.cvboost <- function(object,
                            newx=NULL,
                            ...) {
  if (is.null(newx)) {
    return(object$best_xgb_cvfit$pred)
  }
  else{
    dtest <- xgboost::xgb.DMatrix(data=newx)
    return(predict(object$xgb_fit, newdata=dtest))
  }
}
