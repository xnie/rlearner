#' @include utils.R

# https://topepo.github.io/caret/model-training-and-tuning.html#alternate-performance-metrics
summary_metrics = function(data, lev=NULL, model=NULL) {
	if (is.null(data$weights)) {
		data %<>% dplyr::mutate(weights=1)
	}
	data %>% dplyr::select(obs, pred, weights) %->% c(y, yhat, weights) # makes sure the order is correct
	if (!is.factor(y) && is.numeric(y)) {
		c(wRMSE = sqrt(sum(weights*(y-yhat)^2)/sum(weights)))
	} else {
		positive_class = levels(y)[1]
		y = y==positive_class # convert y to a logical
		p = data[[positive_class]] # the column of probabilities will be called whatever is in levels(y)[1] since y is a factor 
		c(wDeviance = -sum(weights*(y*log(p) + (1-y)*log(1-p)))/sum(weights)) # from the weighted bernoulli log-likelihood
	}
}

# takes a list of fit caret models (hyperparams already optimized by caret), returns the one with lowest weighted MSE (regression) 
# or highest weighted accuracy (classification)
pick_model = function(models) {
	if(length(models)==1) {
		return(models[[1]])
	} else {
		best_model_name = caret::resamples(models)$values %>% # each row is a fold, columns are (model x metric)
		    tidyr::gather(model_metric, value, -Resample) %>% 
		    tidyr::separate(model_metric, c("model","metric"), sep="~") %>%
		    dplyr::group_by(model) %>%
		    dplyr::summarize(mean_value = mean(as.numeric(value), na.rm=T)) %>% # as.numeric in case of weird things because of NAs
		    dplyr::filter(mean_value==ifelse(models[[1]]$maximize, max(mean_value), min(mean_value))) %>%
		    dplyr::pull(model) %>% dplyr::first() # in case of ties
		return(models[[best_model_name]])
	}
}

#' @title Cross-validated supervised learning
#'
#' @details
#' This is a wrapper around \pkg{caret}'s train function, simplifying the interface
#' for the purpose of individual treatment effect meta-learners. Applies cross-validation to select
#' the optimal learning algorithm and hyperparameters. Uses RMSE to select among regression models and 
#' deviance (Bernouilli log-likelihood) to select among probabilistic classifiers.
#' @param x a numeric matrix of features
#' @param y a two-class factor vector for probabilistic classificaton or numeric vector of targets for regression.
#' If y is a factor, the first factor level is treated as the positive class \eqn{c} such that the predicted probabilities 
#' are \eqn{P(Y=c|X)}.
#' @param model_specs a data structure specifying which learning algorithms, hyperparameters should 
#' be cross validated over, and which additional arguments should be 
#' passed to each learner. This should be a list where the names of each element are valid \pkg{caret}
#' methods (learning algorithms). The list element corresponding to each learning algorithm should itself
#' be a list of two elements named \code{tune_grid} (hyperparameters) and \code{extra_args}. 
#' \code{tune_grid} should be a valid \pkg{caret} tune grid of hyperparameters corresponing to the learning algorithm.
#' \code{extra_args} is a named list of additional arguments to be passed on to the learning algorithm. See example.
#' @param weights optional case weights
#' @param k_folds number of cross-validation folds
#' @param select_by optimization method to use for cross-validation: either \code{"best"} for minimum cross-validation
#' error or \code{"oneSE"} for the one-standard-error (1-SE) rule. The implementaion of the 1-SE rule for learners with
#' multiple hyperparameters is governed by \pkg{caret} and may be ad-hoc for some learners. See: \code{\link[caret]{?caret::oneSE}}.
#' @param p_min If provided, probabilities at prediction time will be trimmed to have minimum \code{p_min}. Used for the X-learner.
#' @param p_max If provided, probabilities at prediction time will be trimmed to have maximum \code{p_max}. Used for the X-learner.
#' @examples
#' \dontrun{
#' model_specs = list(
#' gbm = list(
#'     tune_grid = expand.grid(
#'         n.trees = seq(1,501,20), 
#'         interaction.depth=3, 
#'         shrinkage = 0.1, 
#'         n.minobsinnode=3),
#'     extra_args = list(
#'         verbose=F, 
#'         bag.fraction=1)),
#' glmnet = list(
#'     tune_grid = expand.grid(
#'        alpha=c(0,0.5,1),
#'        lambda=exp(seq(-5,2,0.2))),
#'     extra_args = list())
#' )
#' c(x, w, y, ...) %<-% toy_data_simulation(500) # draw a sample 
#' 
#' best_model_y = learner_cv(x, y, model_specs) 
#' y_hat = predict(best_model_y, x)
#' best_model_w = learner_cv(x, w, model_specs)
#' w_hat_prob = predict(best_model_w, x)
#' }
#' @export
learner_cv = function(x, y, model_specs, weights=NULL, k_folds=5, select_by="best", p_min=0, p_max=1) {
	# assume that any binary data coming in is in factor form. 
	# levels(y)[1] should be the name of the positive class (i.e. "treated", or "had_outcome")
	# for instance, to convert a boolean vector w representing treatement to a factor, use:
	# factor(as.factor(w %>% ifelse("treated", "control")), c("treated", "control"))
	if ((select_by=="oneSE") & (length(model_specs)>1)) {
		rlang::abort("The oneSE rule is only defined when comparing models within a single learning algorithm. It is not always clear how to compare the 'complexity' of the models implicit within two different algorihtms (i.e. LASSO and GBM)")
	}
	model = model_specs %>% purrr::imap(function(settings, method) {
		train_args = list(
			x = x, y = y, weights = weights, 
			metric = "wRMSE", maximize=F, # these will be changed if it is a classification problem
			method = method, tuneGrid = settings$tune_grid,
			trControl = trainControl(
				selectionFunction = select_by,
				summaryFunction=summary_metrics,
				method='cv', number=k_folds,
			  	returnResamp="final", savePredictions="final"))
		if(is.factor(y)) {
			train_args$trControl$classProbs=T
			train_args$metric="wDeviance"} # should be minimized just like wRMSE
		do.call(caret::train, c(train_args, settings$extra_args))
	}) %>% pick_model()

	learner = list(model=model)
	if(is.factor(y)) {
		learner$positive_class = levels(y)[1]
		learner$p_min = p_min
		learner$p_max = p_max
	}
	class(learner) = "learner" 
	return(learner)
}

predict.learner = function(object, newdata) {
	if(object$model$modelType == "Classification") {
		predict(object$model, newdata=newdata, type="prob")[[object$positive_class]] %>%
			trim(object$p_min, object$p_max)
	} else {
		predict(object$model, newdata=newdata) 
	}
}

# get the internal cv predictions (at the optimal hyperparameter value) that were used for hyperparameter selection
resample_predictions = function(learner) {
	predictions = learner$model$pred %>% dplyr::arrange(rowIndex)
	if(learner$model$modelType == "Classification") {
		predictions[[learner$positive_class]]
	} else {
		predictions %>% dplyr::pull(pred)
	}
}

#' @title Cross-validated cross-fitting
#'
#' @details
#' Provides both a "deluxe" version and an "economy" version. The deluxe version 
#' preserves data-splitting independence relations. The "economy version" leaks information from the held-out folds 
#' into the predictions on the held-out folds via the hyperparameter selection: data-splitting independence assumptions 
#' do not hold and theoretical guarentees do not follow, but the models fit more quickly. Internal cross validation is
#' performed via \code{\link{learner_cv}}.
#' @param x a numeric matrix of features
#' @param y a two-class factor vector for probabilistic classificaton or numeric vector of targets for regression.
#' If y is a factor, the first factor level is treated as the positive class \eqn{c} such that the predicted probabilities 
#' are \eqn{P(Y=c|X)}.
#' @param model_specs a data structure specifying which learning algorithms, hyperparameters should 
#' be cross validated over, and which additional arguments should be 
#' passed to each learner. This should be a list where the names of each element are valid \pkg{caret}
#' methods (learning algorithms). The list element corresponding to each learning algorithm should itself
#' be a list of two elements named \code{tune_grid} (hyperparameters) and \code{extra_args}. 
#' \code{tune_grid} should be a valid \pkg{caret} tune grid of hyperparameters corresponing to the learning algorithm.
#' \code{extra_args} is a named list of additional arguments to be passed on to the learning algorithm. See example.
#' @param economy flag that determines if "economy" or "deluxe" cross-validated cross-estimation is performed
#' @param weights optional case weights
#' @param k_folds_cf number of cross-fitting folds. Unecessary if \code{economy=T}.
#' @param k_folds number of cross-validation folds
#' @param select_by optimization method to use for cross-validation: either \code{"best"} for minimum cross-validation
#' error or \code{"oneSE"} for the one-standard-error (1-SE) rule. The implementaion of the 1-SE rule for learners with
#' multiple hyperparameters is governed by \pkg{caret} and may be ad-hoc for some learners. See: \code{\link[caret]{?caret::oneSE}}.
#' @examples
#' \dontrun{
#' model_specs = list(
#' gbm = list(
#'     tune_grid = expand.grid(
#'         n.trees = seq(1,501,20), 
#'         interaction.depth=3, 
#'         shrinkage = 0.1, 
#'         n.minobsinnode=3),
#'     extra_args = list(
#'         verbose=F, 
#'         bag.fraction=1)),
#' glmnet = list(
#'     tune_grid = expand.grid(
#'        alpha=c(0,0.5,1),
#'        lambda=exp(seq(-5,2,0.2))),
#'     extra_args = list())
#' )
#' library(zeallot) # imports the %<-% operator, which is syntactic sugar that performs multiple assignment out of a list
#' c(x, w, y, ...) %<-% toy_data_simulation(500) # draw a sample 
#' 
#' y_hat = xval_xfit(x, y, model_specs) 
#' w_hat_prob = xval_xfit(x, w, model_specs)
#' }
#' @export
xval_xfit = function(x, y, model_specs, economy=T, weights=NULL, k_folds_cf=5, k_folds=5, select_by="best") {
	if (economy) {
		learner_cv(x, y, model_specs, weights=weights, 
			k_folds=k_folds, select_by=select_by) %>% 
			resample_predictions()
	} else {
		caret::createFolds(y, k=k_folds_cf) %>%
		purrr::map(function(test_index) {
			learner_cv(x[-test_index,], y[-test_index], model_specs, weights=weights, 
				k_folds=k_folds, select_by=select_by) %>%
				predict(newdata=x[test_index,]) %>%
				data.frame(cross_estimate = ., index=test_index)
		}) %>% dplyr::bind_rows() %>% dplyr::arrange(index) %>% dplyr::pull(cross_estimate)
	}
}
