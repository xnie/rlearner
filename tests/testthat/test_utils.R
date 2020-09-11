context("Toy simulation")
library(rlearner)

n = 100

sanitization_tests = function(x, w, y){
	input =  sanitize_input(x, w, y)
	x = input$x
	w = input$w
	y = input$y
	expect_equal(is.numeric(x), TRUE)
	expect_equal(is.matrix(x), TRUE)
	expect_equal(length(colnames(x)), ncol(x))
	expect_equal(is.logical(w), TRUE)
	expect_equal(is.numeric(y), TRUE)
}

test_that("input sanitization correctly accepts good data and rejects malformed or wrong-type data", {
	sim = toy_data_simulation(n)
	x = sim$x
	w = sim$w
	y = sim$y
	p = sim$p
	m = sim$m
	mu0 = sim$mu0
	mu1 = sim$mu1
	tau = sim$tau


	sanitization_tests(x, w, y)
	sanitization_tests(x, ifelse(w, 1, 0), y)
	x_no_name = x
	colnames(x_no_name) = NULL
	sanitization_tests(x_no_name, w, y)

	expect_error(sanitize_input(x, as.factor(ifelse(w, "t", "c")), y))
	expect_error(sanitize_input(data.frame(x), w, y))
	expect_error(sanitize_input(x, w, y[1:(n-1)]))
	expect_error(sanitize_input(x, w[1:(n-1)], y))
	expect_error(sanitize_input(x[1:(n-1),], w, y))
})
