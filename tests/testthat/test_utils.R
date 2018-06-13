context("Toy simulation")
library(rlearner)
library(zeallot)
library(purrr)

n = 100

sim_tests = function(sim_data) {
	c(x, w, y, p, m, mu0, mu1, tau) %<-% sim_data
	
	# everything is the right size
	expect_equal(
		all(map_int(sim_data[-1], length)==n), # exclude the covariate matrix
		TRUE)
	expect_equal(nrow(x), n)
	expect_equal(
		all(map_lgl(sim_data[-2], is.numeric)), # exclude the treatment vector
		TRUE)
	expect_equal(is.factor(w), TRUE)
}

test_that("toy simulation returns things that are the right size, shape, and type", {
	sim_tests(toy_data_simulation(n)) 
	sim_tests(easy_toy_data_simulation(n))
})
