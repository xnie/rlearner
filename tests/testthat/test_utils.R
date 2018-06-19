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
	expect_equal(is.logical(w), TRUE)
}

test_that("toy simulation returns things that are the right size, shape, and type", {
	sim_tests(toy_data_simulation(n)) 
	sim_tests(easy_toy_data_simulation(n))
})

sanitization_tests = function(x, w, y){
	c(x, w, y) %<-% sanitize_input(x, w, y)
	expect_equal(is.numeric(x), TRUE)
	expect_equal(is.matrix(x), TRUE)
	expect_equal(length(colnames(x)), ncol(x))
	expect_equal(is.logical(w), TRUE)
	expect_equal(is.numeric(y), TRUE)
}

test_that("input sanitization works", {
	c(x, w, y, p, m, mu0, mu1, tau) %<-% toy_data_simulation(n)

	sanitization_tests(x, w, y)
	sanitization_tests(x, ifelse(w, 1, 0), y)

	expect_error(sanitize_input(x, as.factor(ifelse(w, "t", "c")), y))
	expect_error(sanitize_input(x, ifelse(w, 1, 2), y))
	expect_error(sanitize_input(data.frame(x), w, y))
	expect_error(sanitize_input(x, w, w))
	expect_error(sanitize_input(x, w, y[1:(n-1)]))
	expect_error(sanitize_input(x, w[1:(n-1)], y))
	expect_error(sanitize_input(x[1:(n-1),], w, y))
})