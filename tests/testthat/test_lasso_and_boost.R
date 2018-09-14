context("Lasso- and boosting-based learners")
library(rlearner)
library(caret)
library(zeallot)
library(purrr)
library(magrittr)

n = 100
max_tol = 1
min_tol = 1

meta_learner_tests = function(tau_hat, sim_data, mse=0.01) {
  expect_equal(length(tau_hat), length(sim_data$tau))
  expect_equal(is.numeric(tau_hat), TRUE)
  learner_mse = mean((tau_hat - sim_data$tau)^2)
  print(learner_mse)
  expect_equal(learner_mse<mse, TRUE)
}

test_that("lasso based learners return the correct output format and predict well when the problem is easy", {
  set.seed(1)
  easy_sim_data = easy_toy_data_simulation(5*n) # draw a sample
  list(
    rlasso,
    slasso,
    tlasso,
    ulasso,
    xlasso) %>%
    map(function(meta_learner) {
      meta_learner(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y) %>%
        predict(easy_sim_data$x) %>%
        meta_learner_tests(easy_sim_data)
    })
  rs.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=TRUE)
  rs.pred = predict(rs.fit)
  meta_learner_tests(rs.pred, easy_sim_data)
})

test_that("rlasso learns reasonable m_hat on simple data", {
  set.seed(1)
  easy_sim_data = easy_toy_data_simulation(5*n) # draw a sample
  rlasso_fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  print(max(rlasso_fit$m_hat))
  print(max(easy_sim_data$y))
  print(min(rlasso_fit$m_hat))
  print(min(easy_sim_data$y))
  expect_equal((abs(max(rlasso_fit$m_hat) - max(easy_sim_data$y)))/ max(easy_sim_data$y) < max_tol, TRUE)
  expect_equal((abs(min(rlasso_fit$m_hat) - min(easy_sim_data$y)))/ min(easy_sim_data$y) < min_tol, TRUE)
})

test_that("r,s,u,x boosting-based learners product the correct output format and predict well when the problem is easy", {
  set.seed(1)
  easy_sim_data = easy_toy_data_simulation(100*n) # draw a sample
  list(
    rboost,
    sboost,
    uboost,
    xboost) %>%
    map(function(meta_learner) {
      meta_learner(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y) %>%
        predict(easy_sim_data$x) %>%
        meta_learner_tests(easy_sim_data, mse=0.01)
    })
})

test_that("t boosting-based learners product the correct output format and predict well when the problem is easy", {
  set.seed(1)
  easy_sim_data = t_toy_data_simulation(100*n) # draw a sample
  list(tboost) %>%
    map(function(meta_learner) {
      meta_learner(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y) %>%
        predict(easy_sim_data$x) %>%
        meta_learner_tests(easy_sim_data, mse=0.01)
    })
})
