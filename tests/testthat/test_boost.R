context("boosting-based learners")
library(rlearner)

n = 100
max_tol = 1
min_tol = 1

test_that("boosting based learners return the correct output format and predict well when the problem is easy", {
  set.seed(1)
  easy_sim_data = easy_toy_data_simulation(10*n) # draw a sample
  r.fit = rboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  r.pred = predict(r.fit)
  meta_learner_tests(r.pred, easy_sim_data, 0.1)

  r1.fit = rboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1)
  r1.pred = predict(r1.fit)
  invariate_add_tests(r.pred, r1.pred)

  r2.fit = rboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2)
  r2.pred = predict(r2.fit)
  invariate_mult_tests(r.pred, r2.pred)

  s.fit = sboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  s.pred = predict(s.fit)
  meta_learner_tests(s.pred, easy_sim_data, 0.1)

  s1.fit = sboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1)
  s1.pred = predict(s1.fit)
  invariate_add_tests(s.pred, s1.pred)

  s2.fit = sboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2)
  s2.pred = predict(s2.fit)
  invariate_mult_tests(s.pred, s2.pred)

  t.fit = tboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  t.pred = predict(t.fit)
  meta_learner_tests(t.pred, easy_sim_data, 0.1)

  t1.fit = tboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1)
  t1.pred = predict(t1.fit)
  invariate_add_tests(t.pred, t1.pred)

  t2.fit = tboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2)
  t2.pred = predict(t2.fit)
  invariate_mult_tests(t.pred, t2.pred)

  u.fit = uboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  u.pred = predict(u.fit)
  meta_learner_tests(u.pred, easy_sim_data, 0.1)

  u1.fit = uboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1)
  u1.pred = predict(u1.fit)
  invariate_add_tests(u.pred, u1.pred)

  u2.fit = uboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2)
  u2.pred = predict(u2.fit)
  invariate_mult_tests(u.pred, u2.pred)

  x.fit = xboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  x.pred = predict(x.fit)
  meta_learner_tests(x.pred, easy_sim_data, 0.1)

  x1.fit = xboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1)
  x1.pred = predict(x1.fit)
  invariate_add_tests(x.pred, x1.pred)

  x2.fit = xboost(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2)
  x2.pred = predict(x2.fit)
  invariate_mult_tests(x.pred, x2.pred)

})


