context("Lasso-based learners")
library(rlearner)
library(caret)
library(zeallot)
library(purrr)
library(magrittr)

n = 100
max_tol = 1
min_tol = 1

test_that("lasso based r- and s-learners return the correct output format and predict well when the problem has continuous treatments", {
   set.seed(1)
   easy_sim_data = continuous_toy_data_simulation(5*n) # draw a sample

   r.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=F)
   r.pred = predict(r.fit)
   meta_learner_tests(r.pred, easy_sim_data)

   s.fit = slasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
   s.pred = predict(s.fit)
   meta_learner_tests(s.pred, easy_sim_data)

   rs.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=TRUE)
   rs.pred = predict(rs.fit)
   meta_learner_tests(rs.pred, easy_sim_data)
 })

test_that("lasso based learners return the correct output format and predict well when the problem is easy", {
  set.seed(1)
  easy_sim_data = easy_toy_data_simulation(5*n) # draw a sample
  #list(
  #  rlasso,
  #  slasso,
  #  tlasso,
  #  ulasso,
  #  xlasso) %>%
  #  map(function(meta_learner) {
  #    meta_learner(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y) %>%
  #      predict(easy_sim_data$x) %>%
  #      meta_learner_tests(easy_sim_data)
  #  })
  r.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=F)
  r.pred = predict(r.fit)
  meta_learner_tests(r.pred, easy_sim_data)

  s.fit = slasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  s.pred = predict(s.fit)
  meta_learner_tests(s.pred, easy_sim_data)

  t.fit = tlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  t.pred = predict(t.fit)
  meta_learner_tests(t.pred, easy_sim_data)

  u.fit = ulasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  u.pred = predict(u.fit)
  meta_learner_tests(u.pred, easy_sim_data)

  x.fit = xlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  x.pred = predict(x.fit)
  meta_learner_tests(x.pred, easy_sim_data)

  rs.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=TRUE)
  rs.pred = predict(rs.fit)
  meta_learner_tests(rs.pred, easy_sim_data)
})


