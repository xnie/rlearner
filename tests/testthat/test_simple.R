context("simple tests to make sure code runs")
library(rlearner)

n = 50
max_tol = 1
min_tol = 1

test_that("lasso based tests with continuous treatments return the correct output format", {
   set.seed(1)
   easy_sim_data = continuous_toy_data_simulation(n) # draw a sample

   r.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=F)
   r.pred = predict(r.fit)
   simple_meta_learner_tests(r.pred, easy_sim_data)

   s.fit = slasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
   s.pred = predict(s.fit)
   simple_meta_learner_tests(s.pred, easy_sim_data)

   rs.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=TRUE)
   rs.pred = predict(rs.fit)
   simple_meta_learner_tests(rs.pred, easy_sim_data)
 })

test_that("lasso based learners return the correct output format", {
  set.seed(1)
  easy_sim_data = easy_toy_data_simulation(n) # draw a sample
  r.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=F)
  r.pred = predict(r.fit)
  simple_meta_learner_tests(r.pred, easy_sim_data)

  s.fit = slasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  s.pred = predict(s.fit)
  simple_meta_learner_tests(s.pred, easy_sim_data)

  t.fit = tlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  t.pred = predict(t.fit)
  simple_meta_learner_tests(t.pred, easy_sim_data)

  u.fit = ulasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  u.pred = predict(u.fit)
  simple_meta_learner_tests(u.pred, easy_sim_data)

  x.fit = xlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  x.pred = predict(x.fit)
  simple_meta_learner_tests(x.pred, easy_sim_data)

  rs.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=TRUE)
  rs.pred = predict(rs.fit)
  simple_meta_learner_tests(rs.pred, easy_sim_data)
})


test_that("boost based learners return the correct output format", {
  set.seed(1)
  easy_sim_data = easy_toy_data_simulation(n) # draw a sample
  r.fit = rboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  r.pred = predict(r.fit)
  simple_meta_learner_tests(r.pred, easy_sim_data)

  s.fit = sboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  s.pred = predict(s.fit)
  simple_meta_learner_tests(s.pred, easy_sim_data)

  t.fit = tboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  t.pred = predict(t.fit)
  simple_meta_learner_tests(t.pred, easy_sim_data)

  u.fit = uboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  u.pred = predict(u.fit)
  simple_meta_learner_tests(u.pred, easy_sim_data)

  x.fit = xboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  x.pred = predict(x.fit)
  simple_meta_learner_tests(x.pred, easy_sim_data)

})

