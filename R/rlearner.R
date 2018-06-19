#' @include learner_utils.R utils.R

#' @title R-learning for heterogenous treatment effects
#'
#' @details The R-learner estimates heterogenous treatment effects by learning a custom objective function and minimizing it using
#' any suitable machine learning algorithm.
#' @param x a numeric matrix of \strong{covariates}
#' @param w a logical vector indicating \strong{treatment}
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
#' @param k_folds_cf number of cross-estimation folds to use in estimating \code{p_hat} and \code{m_hat}.
#' Unecessary if \code{economy=T} or if both \code{p_hat} and \code{m_hat} are provided.
#' @param k_folds number of cross-validation folds to use in hyperparameter optimization for each model.
#' @param select_by optimization method to use for cross-validation in each model: either \code{"best"} for minimum cross-validation
#' error or \code{"oneSE"} for the one-standard-error (1-SE) rule. The implementaion of the 1-SE rule for learners with
#' multiple hyperparameters is governed by \pkg{caret} and may be ad-hoc for some learners. See: \code{\link[caret]{?caret::oneSE}}.
#' @param p_min If provided, estimated propensities will be trimmed to have minimum \code{p_min}.
#' @param p_max If provided, estimated propensities will be trimmed to have maximum \code{p_max}.
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
#' tau_hat_model = rlearner_cv(x, w, y, model_specs)
#' tau_hat = predict(tau_hat_model, x)
#' }
#' @export
rlearner_cv = function(x, w, y, tau_model_specs,
	p_model_specs=tau_model_specs, m_model_specs=tau_model_specs,
	p_hat=NULL, m_hat=NULL,
	k_folds=5, k_folds_cf=5,
	economy=T, select_by="best",
	p_min=0, p_max=1) {

	c(x, w, y) %<-% sanitize_input(x,w,y)

	if (is.null(p_hat)) {
		p_hat = xval_xfit(x, w, p_model_specs,
			k_folds_cf=k_folds_cf, k_folds=k_folds, economy=economy, select_by=select_by) %>%
			trim(p_min, p_max)
	}
	if (is.null(m_hat)) {
		m_hat = xval_xfit(x, y, m_model_specs,
			k_folds_cf=k_folds_cf, k_folds=k_folds, economy=economy, select_by=select_by)
	}

	r_pseudo_outcome = (y - m_hat)/(w - p_hat)

	r_weights = (w - p_hat)^2

	rlearner = list(
		model=learner_cv(x, r_pseudo_outcome, tau_model_specs, weights=r_weights,
			k_folds=k_folds, select_by=select_by),
		m_hat=m_hat,
		p_hat=p_hat
		)
	class(rlearner) = "rlearner"
	return(rlearner)
}

#' @title Prediction for R-learner
#' @param object a R-learner object
#' @param newx a matrix of covariates for which to predict the treatment effect
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
#' tau_hat_model = rlearner_cv(x, w, y, model_specs)
#' tau_hat = predict(tau_hat_model, x)
#' }
#' @export predict.rlearner
predict.rlearner = function(object, newx, ...) {
	predict(object$model, newdata=newx)
}
