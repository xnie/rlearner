context("Lasso-based learners")
library(rlearner)

n = 100
max_tol = 1
min_tol = 1

test_that("lasso based r- and s-learners return the correct output format and predict well when the problem has continuous treatments", {
   skip_on_cran()
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
  skip_on_cran()
  set.seed(1)
  easy_sim_data = easy_toy_data_simulation(5*n) # draw a sample
  r.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=F)
  r.pred = predict(r.fit)
  meta_learner_tests(r.pred, easy_sim_data)

  r1.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1, rs=F)
  r1.pred = predict(r1.fit)
  invariate_add_tests(r.pred, r1.pred)

  r2.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2, rs=F)
  r2.pred = predict(r2.fit)
  invariate_mult_tests(r.pred, r2.pred)

  s.fit = slasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  s.pred = predict(s.fit)
  meta_learner_tests(s.pred, easy_sim_data)

  s1.fit = slasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1)
  s1.pred = predict(s1.fit)
  invariate_add_tests(s.pred, s1.pred)

  s2.fit = slasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2)
  s2.pred = predict(s2.fit)
  invariate_mult_tests(s.pred, s2.pred)

  t.fit = tlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  t.pred = predict(t.fit)
  meta_learner_tests(t.pred, easy_sim_data)

  t1.fit = tlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1)
  t1.pred = predict(t1.fit)
  invariate_add_tests(t.pred, t1.pred)

  t2.fit = tlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2)
  t2.pred = predict(t2.fit)
  invariate_mult_tests(t.pred, t2.pred)

  u.fit = ulasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  u.pred = predict(u.fit)
  meta_learner_tests(u.pred, easy_sim_data)

  u1.fit = ulasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1)
  u1.pred = predict(u1.fit)
  invariate_add_tests(u.pred, u1.pred, 0.05, 0.01)

  u2.fit = ulasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2)
  u2.pred = predict(u2.fit)
  invariate_mult_tests(u.pred, u2.pred)

  x.fit = xlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  x.pred = predict(x.fit)
  meta_learner_tests(x.pred, easy_sim_data)

  x1.fit = xlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1)
  x1.pred = predict(x1.fit)
  invariate_add_tests(x.pred, x1.pred)

  x2.fit = xlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y*2)
  x2.pred = predict(x2.fit)
  invariate_mult_tests(x.pred, x2.pred)

  rs.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=TRUE)
  rs.pred = predict(rs.fit)
  meta_learner_tests(rs.pred, easy_sim_data)

  rs1.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1, rs=TRUE)
  rs1.pred = predict(rs1.fit)
  invariate_add_tests(rs.pred, rs1.pred)

  rs2.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y+1, rs=TRUE)
  rs2.pred = predict(rs2.fit)
  invariate_mult_tests(rs.pred, rs2.pred)
})

test_that("lasso based learners predict well in setup B in the paper", {
  skip_on_cran()
  set.seed(1)
  easy_sim_data = data_simulation(500) # draw a sample
  r.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=F)
  r.pred = predict(r.fit)
  meta_learner_tests(r.pred, easy_sim_data, 0.5)

  s.fit = slasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  s.pred = predict(s.fit)
  meta_learner_tests(s.pred, easy_sim_data, 0.5)

  t.fit = tlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  t.pred = predict(t.fit)
  meta_learner_tests(t.pred, easy_sim_data, 0.5)

  u.fit = ulasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  u.pred = predict(u.fit)
  meta_learner_tests(u.pred, easy_sim_data, 0.5)

  x.fit = xlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y)
  x.pred = predict(x.fit)
  meta_learner_tests(x.pred, easy_sim_data, 0.5)

  rs.fit = rlasso(easy_sim_data$x, easy_sim_data$w, easy_sim_data$y, rs=TRUE)
  rs.pred = predict(rs.fit)
  meta_learner_tests(rs.pred, easy_sim_data, 0.5)

})



