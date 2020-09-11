context("kernel-based learners")
library(rlearner)

n = 200
max_tol = 1
min_tol = 1

test_that("r,s,u,x kernel-based learners product the correct output format and predict well when the problem is easy", {
  set.seed(1)

 easy_sim_data = data_simulation(n) # draw a sample
  r.fit = rkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  r.pred = predict(r.fit)
  meta_learner_tests(r.pred, easy_sim_data, 0.1)

  s.fit = skern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  s.pred = predict(s.fit)
  meta_learner_tests(s.pred, easy_sim_data, 0.25)
#
  #t.fit = tkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  #t.pred = predict(t.fit)
  #meta_learner_tests(t.pred, easy_sim_data, 0.25)
#
  #u.fit = ukern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  #u.pred = predict(u.fit)
  #meta_learner_tests(u.pred, easy_sim_data, 0.15)
#
  #x.fit = xkern(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  #x.pred = predict(x.fit)
  #meta_learner_tests(x.pred, easy_sim_data, 0.15)

})

