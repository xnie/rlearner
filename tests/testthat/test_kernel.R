context("kernel-based learners")
library(rlearner)
library(caret)
library(zeallot)
library(purrr)
library(magrittr)

n = 100
max_tol = 1
min_tol = 1

test_that("r,s,u,x kernel-based learners product the correct output format and predict well when the problem is easy", {
  set.seed(1)
  easy_sim_data = easy_toy_data_simulation(5*n) # draw a sample
  list(
    rkern,
    skern,
    ukern,
    tkern,
    xkern) %>%
    map(function(meta_learner) {
      meta_learner(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y) %>%
        predict(easy_sim_data$x) %>%
        meta_learner_tests(easy_sim_data, mse=0.1)
    })
})
