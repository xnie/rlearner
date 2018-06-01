# No scaling or centering needs to be done for the U-learner- 
# let each model take care of that internally as per its own defaults

#' @export
U_learner_cv = function(x, w, y, tau_model_specs,
	p_model_specs=NULL, m_model_specs=NULL, 
	p_hat=NULL, m_hat=NULL, 
	k_folds_cv=5, k_folds_ce=5, 
	economy=T, select_by="best",
	p_min=0, p_max=1) {

	if (is.null(p_hat)) {	
		p_hat = xval_xfit(x, w, p_model_specs, 
			k_folds_ce=k_folds_ce, k_folds_cv=k_folds_cv, economy=economy, select_by=select_by) %>%
			trim(p_min, p_max)
	} 
	if (is.null(m_hat)) {
		m_hat = xval_xfit(x, y, mu_model_specs, 
			k_folds_ce=k_folds_ce, k_folds_cv=k_folds_cv, economy=economy, select_by=select_by)
	}

	if (is.factor(w)) {w = w==levels(w)[1]} # turn factor to a logical (the first factor level should be the "treated")
	r_pseudo_outcome = (y - m_hat)/(w - p_hat)

	U_learner = list(
		model=learner_cv(x, r_pseudo_outcome, tau_model_specs, 
			k_folds=k_folds_cv, select_by=select_by) 
		)
	class(U_learner) = "U_learner"
	return(U_learner)
}

#' @export
predict.U_learner = function(object, x) {
	predict(object$model, newdata=x)
}