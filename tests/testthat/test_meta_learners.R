context("Meta-learners")
library(rlearner)
library(caret)
library(zeallot)
library(purrr)
library(magrittr)

n = 100
#  # draw a sample 
# weights = runif(n)

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

meta_learner_tests = function(tau_hat, sim_data) {
	expect_equal(length(tau_hat), length(sim_data$tau))
	expect_equal(is.numeric(tau_hat), TRUE)
	mse = mean((tau_hat - sim_data$tau)^2)
	expect_equal(mse<0.01, TRUE)
}

test_that("meta-learners return sensible things and predict well when the problem is easy", {
	set.seed(1)
	easy_sim_data = easy_toy_data_simulation(5*n) # draw a sample 
	list(
		rlearner_cv, 
		slearner_cv, 
		tlearner_cv,
		ulearner_cv, 
		xlearner_cv) %>%
	map(function(meta_learner) {
		meta_learner(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, model_specs) %>% 
		predict(easy_sim_data$x) %>%
		meta_learner_tests(easy_sim_data)
	})
})

test_that("R-learner internals are good", {
	c(x, w, y, p, m, mu0, mu1, tau) %<-% toy_data_simulation(n)
	model = rlearner_cv(x, w, y, model_specs)
	expect_equal(class(model), "rlearner")
	expect_equal(class(model$model), "learner")
	expect_equal(length(model$m_hat), nrow(x))
	expect_equal(length(model$p_hat), nrow(x))
	expect_equal(is.numeric(model$m_hat), TRUE)
	expect_equal(is.numeric(model$p_hat), TRUE)
})

test_that("S-learner internals are good", {
	c(x, w, y, p, m, mu0, mu1, tau) %<-% toy_data_simulation(n)
	model = slearner_cv(x, w, y, model_specs)
	expect_equal(class(model), "slearner")
	expect_equal(class(model$model), "learner")
	expect_equal(class(model$standardization), "preProcess")
})

test_that("T-learner internals are good", {
	c(x, w, y, p, m, mu0, mu1, tau) %<-% toy_data_simulation(n)
	model = tlearner_cv(x, w, y, model_specs)
	expect_equal(class(model), "tlearner")
	model %>% map(~expect_equal(class(.), "learner"))
})

test_that("U-learner internals are good", {
	c(x, w, y, p, m, mu0, mu1, tau) %<-% toy_data_simulation(n)
	model = ulearner_cv(x, w, y, model_specs)
	expect_equal(class(model), "ulearner")
	expect_equal(class(model$model), "learner")
})

test_that("X-learner internals are good", {
	c(x, w, y, p, m, mu0, mu1, tau) %<-% toy_data_simulation(n)
	model = xlearner_cv(x, w, y, model_specs)
	expect_equal(class(model), "xlearner")
	expect_equal(class(model$p_hat_model), "learner")
	expect_equal(class(model$tau0_hat_model), "learner")
	expect_equal(class(model$tau1_hat_model), "learner")
})

