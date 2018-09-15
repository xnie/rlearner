#' @include learner_utils.R utils.R
#'
#' @title X-learner using generic black-box machine learning model from caret
#'
#' @description X-learner as proposed by Künzel, Sekhon, Bickel, and Yu (2017), using generic black-box machine learning model from caret
#'
#' @param x a numeric matrix of \strong{covariates}
#' @param w a logical vector indicating \strong{treatment}
#' @param y a numeric vector of \strong{outcomes}
#' @param tau_model_specs specification for the model of \eqn{\tau(x) = E[Y(1) - Y(0)|X=x]}. See \code{\link{learner_cv}}.
#' @param p_model_specs specification for the model of \eqn{p(x) = E[W|X=x]}. See \code{\link{learner_cv}}.
#' Not needed if \code{p_hat} is provided.
#' @param mu0_hat a numeric vector of estimates of the outcome under control for all observations.
#' The X-learner will estimate these values internally if not provided.
#' @param mu1_hat a numeric vector of estimates of the outcome under treatment for all observations.
#' The X-learner will estimate these values internally if not provided.
#' @param mu_model_specs specification for the models of \eqn{m_w(x) = E[Y|W=w,X=x]}. See \code{\link{learner_cv}}.
#' Not needed if \code{mu0_hat_0} and \code{mu0_hat} are provided.
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
#' tau_hat_model = xlearner_cv(x, w, y, model_specs)
#' tau_hat = predict(tau_hat_model, x)
#' }
#' @export
xlearner_cv = function(x, w, y, tau_model_specs,
	p_model_specs=tau_model_specs, mu_model_specs=tau_model_specs,
	mu0_hat=NULL, mu1_hat=NULL,
	k_folds=5, select_by="best",
	p_min=0, p_max=1) {

	c(x, w, y) %<-% sanitize_input(x,w,y)

	p_hat_model = learner_cv(x, w, p_model_specs,
		k_folds=k_folds, select_by=select_by,
		p_min = p_min, p_max=p_max)

	if (is.null(mu1_hat)) {
		mu1_hat_model = learner_cv(x[w,], y[w], mu_model_specs,
			k_folds=k_folds, select_by=select_by)
		mu1_hat = predict(mu1_hat_model, x)
	}

	if (is.null(mu0_hat)) {
		mu0_hat_model = learner_cv(x[!w,], y[!w], mu_model_specs,
			k_folds=k_folds, select_by=select_by)
		mu0_hat = predict(mu0_hat_model, x)
	}

	d1 = y[w] - mu0_hat[w]
	d0 = mu1_hat[!w] - y[!w]

	tau1_hat_model = learner_cv(x[w,], d1, tau_model_specs,
			k_folds=k_folds, select_by=select_by)
	tau0_hat_model = learner_cv(x[!w,], d0, tau_model_specs,
		k_folds=k_folds, select_by=select_by)

	xlearner = list(
		p_hat_model = p_hat_model,
		tau0_hat_model = tau0_hat_model,
		tau1_hat_model = tau1_hat_model)
	class(xlearner) = "xlearner"
	return(xlearner)
}

#' @title Prediction for X-learner
#' @param object a X-learner object
#' @param newx a matrix of covariates for which to predict the treatment effect
#' @param ... additional arguments (currently not used)
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
#' tau_hat_model = xlearner_cv(x, w, y, model_specs)
#' tau_hat = predict(tau_hat_model, x)
#' }
#' @export predict.xlearner
predict.xlearner = function(object, newx, ...) {
  newx = sanitize_x(newx)
	object %>% purrr::map(function(model) {
		predict(model, newx=newx)
	}) %->%	c(p_hat, tau0_hat, tau1_hat)
	p_hat = trim(p_hat, object$p_hat_model$p_min, object$p_hat_model$p_max)
	return(p_hat*tau0_hat + (1-p_hat)*tau1_hat)
}
