

X_learner_cv = function(x, w, y, tau_model_specs,
	p_model_specs, mu_model_specs=NULL, 
	mu0_hat_1=NULL, mu1_hat_0=NULL,
	k_folds=5, select_by="best") {
	
	p_hat_model = learner_cv(x, w, p_model_specs, 
		k_folds=k_folds, select_by=select_by)

	w = w==levels(w)[1]

	if (is.null(mu1_hat_0)) {
		mu1_hat_model = learner_cv(x[w,], y[w], mu_model_specs, 
			k_folds=k_folds, select_by=select_by)
		mu1_hat_0 = predict(mu1_hat_model, x[!w,])
	}

	if (is.null(mu0_hat_1)) {
		mu0_hat_model = learner_cv(x[!w,], y[!w], mu_model_specs, 
			k_folds=k_folds, select_by=select_by)
		mu0_hat_1 = predict(mu0_hat_model, x[w,])
	}

	d1 = y[w] - mu0_hat_1
	d0 = mu1_hat_0 - y[!w]

	tau1_hat_model = learner_cv(x[w,], d1, tau_model_specs, 
			k_folds=k_folds, select_by=select_by)
	tau0_hat_model = learner_cv(x[!w,], d0, tau_model_specs, 
		k_folds=k_folds, select_by=select_by)

	X_learner = list(
		p_hat_model = p_hat_model,
		tau0_hat_model = tau0_hat_model,
		tau1_hat_model = tau1_hat_model)
	class(X_learner) = "X_learner"
	return(X_learner)
}

predict.X_learner = function(object, x, p_min=0, p_max=1) {
	object %>% map(function(model) {
		predict(model, x)
	}) %->%	c(p_hat, tau0_hat, tau1_hat)
	p_hat = threshold(p_hat, p_min, p_max)
	return(p_hat*tau0_hat + (1-p_hat)*tau1_hat)
}