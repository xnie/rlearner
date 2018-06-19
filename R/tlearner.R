#' @include learner_utils.R utils.R

#' @title T-learning for heterogenous treatment effects
#'
#' @param x a numeric matrix of \strong{covariates}
#' @param w a logical vector indicating \strong{treatment}
#' @param y a numeric vector of \strong{outcomes}
#' @param model_specs specification for the models of \eqn{\mu_w(x) = E[Y|W=w,X=x]}. See \code{\link{learner_cv}}.
#' @param k_folds number of cross-validation folds to use in hyperparameter optimization for each model.
#' @param select_by optimization method to use for cross-validation in each model: either \code{"best"} for minimum cross-validation
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
#' tau_hat_model = tlearner_cv(x, w, y, model_specs)
#' tau_hat = predict(tau_hat_model, x)
#' }
#' @export
tlearner_cv = function(x, w, y, model_specs, k_folds=5, select_by="best") {
	c(x, w, y) %<-% sanitize_input(x,w,y)

	tlearner = c(TRUE, FALSE) %>% map(function(condition) {
		learner_cv(
			x[w==condition,], y[w==condition], # filters to only include subjects treated under "condition"
			model_specs, k_folds=k_folds)
	}) # this is a list of two fit models
	class(tlearner) = "tlearner"
	return(tlearner)
}

#' @title Prediction for U-learner
#' @param object a U-learner object
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
#' tau_hat_model = tlearner_cv(x, w, y, model_specs)
#' tau_hat = predict(tau_hat_model, x)
#' }
#' @export predict.tlearner
predict.tlearner = function(object, newx, ...) {
	object %>%
		map(~predict(., newdata=newx)) %->%
		c(mu1_hat, mu0_hat) # these will come out in this order because of the above: c(T,F) %>% map(...
	return(mu1_hat - mu0_hat)
}
