context("kernel-based learners")
library(rlearner)

n = 500
max_tol = 1
min_tol = 1

test_that("r,s,u,x kernel-based learners produce expected performance in setup B in the paper", {
  set.seed(1)

  easy_sim_data = data_simulation(n) # draw a sample
  r.fit = rkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y,
                         k_folds = 5)
  r.pred = predict(r.fit)
  meta_learner_tests(r.pred, easy_sim_data, 0.1)

  r1.fit = rkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1,
                         k_folds = 5)
  r1.pred = predict(r1.fit)
  invariate_add_tests(r.pred, r1.pred)

  r2.fit = rkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2,
                         k_folds = 5)
  r2.pred = predict(r2.fit)
  invariate_mult_tests(r.pred, r2.pred)

  s.fit = skern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, k_folds=5)
  s.pred = predict(s.fit)
  meta_learner_tests(s.pred, easy_sim_data, 0.15)

  s1.fit = skern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1, k_folds=5)
  s1.pred = predict(s1.fit)
  invariate_add_tests(s.pred, s1.pred)

  s2.fit = skern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2, k_folds=5)
  s2.pred = predict(s2.fit)
  invariate_mult_tests(s.pred, s2.pred)

  t.fit = tkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, k_folds=5)
  t.pred = predict(t.fit)
  meta_learner_tests(t.pred, easy_sim_data, 0.15)

  t1.fit = tkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1, k_folds=5)
  t1.pred = predict(t1.fit)
  invariate_add_tests(t.pred, t1.pred)

  t2.fit = tkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2, k_folds=5)
  t2.pred = predict(t2.fit)
  invariate_mult_tests(t.pred, t2.pred)

  u.fit = ukern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, k_folds=5)
  u.pred = predict(u.fit)
  meta_learner_tests(u.pred, easy_sim_data, 0.1)

  u1.fit = ukern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1, k_folds=5)
  u1.pred = predict(u1.fit)
  invariate_add_tests(u.pred, u1.pred)

  u2.fit = ukern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2, k_folds=5)
  u2.pred = predict(u2.fit)
  invariate_mult_tests(u.pred, u2.pred)

  x.fit = xkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, k_folds=5)
  x.pred = predict(x.fit)
  meta_learner_tests(x.pred, easy_sim_data, 0.1)

  x1.fit = xkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1, k_folds=5)
  x1.pred = predict(x1.fit)
  invariate_add_tests(x.pred, x1.pred)

  x2.fit = xkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2, k_folds=5)
  x2.pred = predict(x2.fit)
  invariate_mult_tests(x.pred, x2.pred)

})

