# sanitize input to xval_xfit

#' @import magrittr
#' @import caret
#' @import zeallot
#' @import tidyverse

# All of these functions are fancy wrappers around calls to caret

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
		p = data %>% pull(!!positive_class) # the column of probabilities will be called whatever is in levels(y)[1] since y is a factor 
		c(wDeviance = -sum(weights*(y*log(p) + (1-y)*log(1-p)))/sum(weights)) # from the weighted bernoulli log-likelihood
	}
}

# takes a list of fit caret models (hyperparams already optimized by caret), returns the one with lowest weighted MSE (regression) 
# or highest weighted accuracy (classification)
pick_model = function(models) {
	if(length(models)==1) {
		return(models[[1]])
	} else {
		best_model_name = resamples(models)$values %>% # each row is a fold, columns are (model x metric)
		    tidyr::gather(model_metric, value, -Resample) %>% 
		    tidyr::separate(model_metric, c("model","metric"), sep="~") %>%
		    dplyr::group_by(model) %>%
		    dplyr::summarize(mean_value = mean(as.numeric(value), na.rm=T)) %>% # as.numeric in case of weird things because of NAs
		    dplyr::filter(mean_value==ifelse(models[[1]]$maximize, max(mean_value), min(mean_value))) %>%
		    dplyr::pull(model) %>% dplyr::first() # in case of ties
		return(models[[best_model_name]])
	}
}

# fits models via caret, uses CV to select one, returns it
# selection should be either "best" or "oneSE"
learner_cv = function(x, y, model_specs, weights=NULL, k_folds=5, select_by="best") {
	# assume that any binary data coming in is in factor form. 
	# levels(y)[1] should be the name of the positive class (i.e. "treated", or "had_outcome")
	# for instance, to convert a boolean vector w representing treatement to a factor, use:
	# factor(as.factor(w %>% ifelse("treated", "control")), c("treated", "control"))
	if(select_by=="oneSE") {
		rlang::warn("The oneSE rule uses heuristics set in the caret package to rank models by complexity. See ?caret::oneSE")
		if(length(model_specs)>1) {
			rlang::abort("The oneSE rule is only defined when comparing models within a single learning algorithm. It is not always clear how to compare the 'complexity' of the models implicit within two different algorihtms (i.e. LASSO and GBM)")}}
	model = model_specs %>% imap(function(settings, method) {
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
		do.call(train, c(train_args, settings$extra_args))
	}) %>% pick_model()

	learner = list(model=model)
	if(is.factor(y)) {
		learner$positive_class = levels(y)[1]
	}
	class(learner) = "learner" 
	return(learner)
}

predict.learner = function(object, newdata) {
	if(object$model$modelType == "Classification") {
		predict(object$model, newdata=newdata, type="prob")[[object$positive_class]]
	} else {
		predict(object$model, newdata=newdata) 
	}
}

# get the internal cv predictions (at the optimal hyperparameter value) that were used for hyperparameter selection
resample_predictions = function(learner) {
	predictions = learner$model$pred %>% arrange(rowIndex)
	if(learner$model$modelType == "Classification") {
		positive_class = learner$positive_class
		predictions %>% pull(!!positive_class)
	} else {
		predictions %>% pull(pred)
	}
}

#' Cross-validated cross-estimation. Provides both a "deluxe" version and an "economy" version. The deluxe version 
#' preserves data-splitting independence relations. The "economy version" leaks information from the held-out folds 
#' into the predictions on the held-out folds via the hyperparameter selection. Data-splitting independence assumptions 
#' do not hold and theoretical guarentees do not follow, but the models fit more quickly
#' @export
xval_xfit = function(x, y, model_specs, economy=T, weights=NULL, k_folds_ce=5, k_folds_cv=5, select_by="best") {
	if (economy) {
		learner_cv(x, y, model_specs, weights=weights, 
			k_folds=k_folds_cv, select_by=select_by) %>% 
			resample_predictions()
	} else {
		createFolds(y, k=k_folds_ce) %>%
		map(function(test_index) {
			learner_cv(x[-test_index,], y[-test_index], model_specs, weights=weights, 
				k_folds=k_folds_cv, select_by=select_by) %>%
				predict(newdata=x[test_index,]) %>%
				data.frame(cross_estimate = ., index=test_index)
		}) %>% bind_rows %>% dplyr::arrange(index) %>% dplyr::pull(cross_estimate)
	}
}
