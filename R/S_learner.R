
#' @export
S_learner_cv = function(x, w, y, model_specs, k_folds=5, select_by="best") {
	if (is.factor(w)) {w = w==levels(w)[1]} # turn factor to a logical (the first factor level should be the "treated")

	if ("glmnet" %in% names(model_specs)) { # tell glmnet not to standardize... other models may also be standardizing so caveat emptor
		model_specs$glmnet$extra_args$standardize = F
	}

	standardization = caret::preProcess(x, method=c("center", "scale")) # get the standardization params
	x = predict(standardization, x)							 # standardize the input
	x_expanded = cbind(x, (w-0.5)*x, (w-0.5)) 
# check that the names don't mess things up
	# it's not clear how, in general, to have different regularization on x and (w-0.5)x, so the "fancy" S-learner
	# is difficult to implement in a general purpose way.
	# note that glmnet will add its own intercept and won't regularize it

	S_learner = list(
		model = learner_cv(x_expanded, y, model_specs, k_folds=k_folds, select_by=select_by),
		standardization = standardization) 
	class(S_learner) = "S_learner"
	return(S_learner)
}

#' @export
predict.S_learner = function(object, x) {
	x = predict(object$standardization, x) # standardize the new data using the same standardization as with the training data
	list(0, 1) %>% purrr::map(function(w) {
		predict(object$model, newdata=cbind(x, (w-0.5)*x, (w-0.5)))
	}) %->% c(mu0_hat, mu1_hat)
	return(mu1_hat - mu0_hat)
}
