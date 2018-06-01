
#' @title R-learning for heterogenous treatment effects
#'
#' @details The R-learner estimates heterogenous treatment effects by learning a custom objective function and minimizing it using
#' any suitable machine learning algorithm. 
#' @param x a numeric matrix of \strong{covariates}
#' @param w a two-class factor vector of \strong{treatments}. The first factor level is treated as the positive class \eqn{w=1}
#' @param y a numeric vector of \strong{outcomes}
#' @param tau_model_specs specification for the model of \eqn{\tau(x) = E[Y(1) - Y(0)|X=x]}. See \code{\link{learner_cv}}.
#' @param p_hat a numeric vector of estimates of the treatment propensity \eqn{p(x) = E[W|X=x]} of each observation. 
#' The R-learner will estimate these values internally using cross-validated cross-estimation if \code{p_hat} is not provided.
#' @param m_hat a numeric vector of estimates of the outcome marginalized over the treatment (\eqn{m(x) = E[Y|X=x]}) for each observation.
#' The R-learner will estimate these values internally using cross-validated cross-estimation if \code{m_hat} is not provided.
#' @param p_model_specs specification for the model of \eqn{p(x) = E[W|X=x]}. See \code{\link{learner_cv}}.
#' Not needed if \code{p_hat} is provided.
#' @param m_model_specs specification for the model of \eqn{m(x) = E[Y|X=x]}. See \code{\link{learner_cv}}.
#' Not needed if \code{m_hat} is provided.
#' @param economy flag that determines if "economy" or "deluxe" cross-validated cross-estimation is performed when learning
#' the models for \eqn{p(x)} and \eqn{m(x)}.
#' Not needed if both \code{p_hat} and \code{m_hat} are provided.
#' @param k_folds_ce number of cross-estimation folds to use in estimating \code{p_hat} and \code{m_hat}.
#' Unecessary if \code{economy=T} or if both \code{p_hat} and \code{m_hat} are provided.
#' @param k_folds_cv number of cross-validation folds to use in hyperparameter optimization for each model.
#' @param select_by optimization method to use for cross-validation in each model: either \code{"best"} for minimum cross-validation
#' error or \code{"oneSE"} for the one-standard-error (1-SE) rule. The implementaion of the 1-SE rule for learners with
#' multiple hyperparameters is governed by \pkg{caret} and may be ad-hoc for some learners. See: \code{\link[caret]{?caret::oneSE}}.
#' @param p_min If provided, estiamted propensities will be trimmed to have minimum \code{p_min}.
#' @param p_max If provided, estiamted propensities will be trimmed to have maximum \code{p_max}.
#' @return 
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
#' n = 500
#' x = data.frame("covariate_1" = rnorm(n), "covariate_2"= rnorm(n)) %>% 
#'     make_matrix
#' logit_p = (x %*% c(1,1))
#' p = exp(logit_p)/(1+exp(logit_p))
#' w = rbinom(n,1,p)==1
#' w_factor = factor(as.factor(w %>% ifelse("treated", "control")), c("treated", "control"))
#' tau = (x %*% c(1,1))^2
#' m = x %*% c(1,-3)
#' y = (m + tau/2*(2*w-1))[,1]
#' 
#' tau_hat_model = R_learner_cv(x, y, model_specs) 
#' tau_hat = predict(tau_hat_model, x)
#' }
#' @export
R_learner_cv = function(x, w, y, tau_model_specs,
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
		m_hat = xval_xfit(x, y, m_model_specs, 
			k_folds_ce=k_folds_ce, k_folds_cv=k_folds_cv, economy=economy, select_by=select_by)
	}

	if (is.factor(w)) {w = w==levels(w)[1]} # turn factor to a logical (the first factor level should be the "treated")
	r_pseudo_outcome = (y - m_hat)/(w - p_hat)
	r_weights = (w - p_hat)^2

	R_learner = list(
		model=learner_cv(x, r_pseudo_outcome, tau_model_specs, weights=r_weights, 
			k_folds=k_folds_cv, select_by=select_by) 
		)
	class(R_learner) = "R_learner"
	return(R_learner)
}

#' @export
predict.R_learner = function(object, x) {
	predict(object$model, newdata=x)
}
