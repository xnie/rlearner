
# No scaling or centering needs to be done for the T-learner- 
# let each model take care of that internally as per its own defaults

#' @export
T_learner_cv = function(x, w, y, model_specs, k_folds=4, select_by="best") {
	T_learner = levels(w) %>% map(function(condition) {
		learner_cv(
			x[w==condition,], y[w==condition], # filters to only include subjects treated under "condition"
			model_specs, k_folds=k_folds)
	}) # this is a list of two fit models
	class(T_learner) = "T_learner"
	return(T_learner)
}

#' @export
predict.T_learner = function(object, x) {
	object %>% 
		map(~predict(., newdata=x)) %->% 
		c(mu1_hat, mu0_hat) # these will come out in the order of levels(w), where the first level indicates positive treatment
	return(mu1_hat - mu0_hat)
}
