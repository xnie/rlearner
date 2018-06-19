context("Interface to caret")
library(rlearner)
library(caret)
library(zeallot)
library(purrr)
library(magrittr)

n = 100
c(x, w, y, p, m, mu0, mu1, tau) %<-% toy_data_simulation(n) # draw a sample
weights = runif(n)

model_specs = list(
	gbm = list(
	    tune_grid = expand.grid(
	        n.trees = seq(1,501,20),
	        interaction.depth=3,
	        shrinkage = 0.1,
	        n.minobsinnode=3),
	    extra_args = list(
	        verbose=F,
	        bag.fraction=1)),
	glmnet = list(
	    tune_grid = expand.grid(
	       alpha=c(0,0.5,1),
	       lambda=exp(seq(-5,2,0.2))),
	    extra_args = list())
)

reg_tests = function(reg_model) {
	expect_equal(class(reg_model), "learner")
	expect_equal(class(reg_model$model), "train")
	expect_equal(reg_model$model$modelType, "Regression")
	expect_equal(reg_model$model$maximize, FALSE)
	expect_equal(reg_model$model$metric, "wRMSE")
	expect_equal(nrow(reg_model$model$pred), n)
	y_hat = predict(reg_model, x)
	expect_equal(length(y_hat), nrow(x))
	expect_equal(is.numeric(y_hat), TRUE)
}

cls_tests = function(cls_model, p_min=0, p_max=1) {
	expect_equal(cls_model$model$modelType, "Classification")
	expect_equal(cls_model$model$maximize, FALSE)
	expect_equal(cls_model$model$metric, "wDeviance")
	expect_equal(all(levels(w) %in% names(cls_model$model$pred)), TRUE)
	expect_equal(nrow(cls_model$model$pred), n)
	w_hat = predict(cls_model, x)
	expect_equal(length(w_hat), nrow(x))
	expect_equal(is.numeric(w_hat), TRUE)
	expect_equal(all((p_min <= w_hat) & (w_hat <= p_max)), TRUE)
}

test_that("caret returns a correctly-structured and functional trained model object when doing cross-validated regression", {
	reg_tests(learner_cv(x, y, model_specs))
	reg_tests(learner_cv(x, y, model_specs, weights=weights))
	expect_error(learner_cv(x, y, model_specs, select_by="oneSE"))
	reg_tests(learner_cv(x, y, model_specs[1], select_by="oneSE"))
})

test_that("caret returns a correctly-structured and functional trained model object when doing cross-validated probabilistic classification", {
	cls_tests(learner_cv(x, w, model_specs))
	cls_tests(learner_cv(x, w, model_specs, weights=weights))
	c(p_min, p_max) %<-% list(0.45, 0.55)
	cls_tests(learner_cv(x, w, model_specs, p_min=p_min, p_max=p_max),
		p_min=p_min, p_max=p_max)
})

test_that("the conditional mean outcome can be easily predicted when the true model is linear, there is no noise, and there are many samples", {
	set.seed(1)
	c(x, w, y, p, m, mu0, mu1, tau) %<-% easy_toy_data_simulation(5*n) # draw a sample
	y_hat = learner_cv(x, y, model_specs) %>% predict(x)
	mse = mean((y_hat - m)^2)
	expect_equal(mse<0.1, TRUE)
})

xval_xfit_tests = function(est) {
	expect_equal(length(est), nrow(x))
	expect_equal(is.numeric(est), TRUE)
}

test_that("cross-validated cross-fitting returns vectors of predictions that are the right length and type", {
	list(
		xval_xfit(x, y, model_specs),
		xval_xfit(x, w, model_specs, economy=F),
		xval_xfit(x, y, model_specs[2], select_by="oneSE")
	) %>% map(xval_xfit_tests)
})

