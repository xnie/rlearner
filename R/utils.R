# For thresholding propensity scores
trim = function(x, min, max) {
	x[x>max] = max
	x[x<min] = min
	return(x)
}

sanitize_x = function(x){
	# make sure x is a numeric matrix with named columns (for caret)
	if (!is.matrix(x) | !is.numeric(x) | any(is.na(x))) {
		stop("x must be a numeric matrix with no missing values")
	}
	colnames(x) = stringr::str_c("covariate_", 1:ncol(x))
	return(x)
}

sanitize_input = function(x,w,y) {
  x = sanitize_x(x)

  if (!is.numeric(w)) {
		stop("the input w should be a numeric vector")
  }
	if (is.numeric(w) & all(w %in% c(0,1))) {
		w = w==1
	}

	# make sure y is a numeric vector
	if (!is.numeric(y)) {
		stop("y should be a numeric vector")
	}

	# make sure the dimensions align
	if (length(y)!=nrow(x) | length(w)!=nrow(x)) {
		stop("nrow(x), length(w), and length(y) should all be equal")
	}

	return(list(x=x,
	            w=w,
	            y=y))
}

#' @title data simulation
#'
#' @description Generates a dataset of size \eqn{n} that can be used to experiment with the learners and meta-learners
#' in this package.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' data = data_simulation(500) # draw a sample
#' @export
data_simulation = function(n) {
  x = stats::model.matrix(~.-1, data.frame("covariate_1" = stats::rnorm(n), "covariate_2"= stats::rnorm(n), "covariate_3" = stats::rnorm(n), "covariate_4" = stats::rnorm(n), "covariate_5" = stats::rnorm(n), "covariate_6" = stats::rnorm(n)))
	p = 0.5
	w = as.numeric(stats::rbinom(n,1,p)==1)
  m = pmax(0, x[,1] + x[,2], x[,3]) + pmax(0, x[,4] + x[,5])
  tau = x[,1] + log(1 + exp(x[,2]))
	mu1 = m + tau/2
	mu0 = m - tau/2
	y = w*mu1 + (1-w) * mu0 + 0.5*stats::rnorm(n)
	list(x=x, w=w, y=y, p=p, m=m, mu0=mu0, mu1=mu1, tau=tau)
}

#' @title Toy data simulation
#'
#' @description Generates a toy dataset of size \eqn{n} that can be used to experiment with the learners and meta-learners
#' in this package. The generative process should be easy to learn with linear methods.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' data = easy_toy_data_simulation(500) # draw a sample
#' @export
easy_toy_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = stats::rnorm(n), "covariate_2"= stats::rnorm(n)))
	p = rep(0.5, n)
	w = as.numeric(stats::rbinom(n,1,p)==1)
	tau = x %*% c(1,1)
	m = x %*% c(0.5,-0.5)
	mu1 = m + tau/2
	mu0 = m - tau/2
	y = (m + tau/2*(2*w-1))[,1]
	list(x=x, w=w, y=y, p=p, m=m, mu0=mu0, mu1=mu1, tau=tau)
}
#' @title Toy data simulation with continuous treatments
#'
#' @description Generates a toy dataset of size \eqn{n} that can be used to experiment with the learners and meta-learners
#' in this package. The generative process should be easy to learn with linear methods.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' data = continuous_toy_data_simulation(500) # draw a sample
#' @export
continuous_toy_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = stats::rnorm(n), "covariate_2"= stats::rnorm(n)))
	w = stats::runif(n, 0, 1)
	tau = x %*% c(1,1)
	m = x %*% c(0.5,-0.5)
	mu1 = m + tau/2
	mu0 = m - tau/2
	y = (m + tau/2*(2*w-1))[,1]
	list(x=x, w=w, y=y,  m=m, mu0=mu0, mu1=mu1, tau=tau)
}

#' @title Toy data simulation for T-learner
#'
#' @description Generates a toy dataset of size \eqn{n} that can be used to experiment with T-learners
#' in this package. The generative process should be easy to learn with linear methods.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' data = t_toy_data_simulation(500) # draw a sample
#' @export
t_toy_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = stats::rnorm(n), "covariate_2"= stats::rnorm(n)))
	p = rep(0.5, n)
	w = as.numeric(stats::rbinom(n,1,p)==1)
	mu1 = sin(x[,1] * 2)
	mu0 = x[,2] * 3 + 10
	y = w * mu1 + (1-w) * mu0
	tau = mu1 - mu0
	m = p * mu1 + (1-p) * mu0
	list(x=x, w=w, y=y, p=p, m=m, mu0=mu0, mu1=mu1, tau=tau)
}
#' @title data simulation for T-learner
#'
#' @description Generates a dataset of size \eqn{n} that can be used to experiment with T-learners
#' in this package. The generative process should be easy to learn with linear methods.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' data = t_data_simulation(500) # draw a sample
#' @export
t_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = stats::rnorm(n), "covariate_2"= stats::rnorm(n), "covariate_3" = stats::rnorm(n), "covariate_4" = stats::rnorm(n), "covariate_5" = stats::rnorm(n), "covariate_6" = stats::rnorm(n)))
  p = 1/(1 + exp(-x[,1]) + exp(-x[,2]))
	w = as.numeric(stats::rbinom(n,1,p)==1)
  b = (pmax(x[,1] + x[,2] + x[,3], 0) + pmax(x[,4] + x[,5], 0)) / 2
  tau = pmax(x[,1] + x[,2] + x[,3], 0) - pmax(x[,4] + x[,5], 0)

  mu1 = b + 0.5 * tau
  mu0 = b - 0.5 * tau

	y = w * mu1 + (1-w) * mu0 + stats::rnorm(n)
	m = p * mu1 + (1-p) * mu0
	list(x=x, w=w, y=y, p=p, m=m, mu0=mu0, mu1=mu1, tau=tau)
}
#' @title helper function for result comparison
#' @description helper function to check if the learned tau_hat is within some bounded error from the groundtruth
#' @param tau_hat user-supplied treatment effect estimate
#' @param sim_data simulated groundtruth data
#' @param mse user-supplied error tolerance
#'
#' @export
meta_learner_tests = function(tau_hat, sim_data, mse=0.01) {
  testthat::expect_equal(length(tau_hat), length(sim_data$tau))
  testthat::expect_equal(is.numeric(tau_hat), TRUE)
  learner_mse = mean((tau_hat - sim_data$tau)^2)
  print(learner_mse)
  testthat::expect_equal(learner_mse<mse, TRUE)
}

#' @title helper function for testing the code runs
#' @description helper function to check if the learned tau_hat is numeric
#' @param tau_hat user-supplied treatment effect estimate
#' @param sim_data simulated groundtruth data
#' @param mse user-supplied error tolerance
#'
#' @export
simple_meta_learner_tests = function(tau_hat, sim_data, mse=0.01) {
  testthat::expect_equal(length(tau_hat), length(sim_data$tau))
  testthat::expect_equal(is.numeric(tau_hat), TRUE)
}

#' @title helper function for testing treatment effect is invariant when outcome adds 1
#' @description helper function to test treatment effect is invariant when outcome adds 1
#' @param tau_hat user-supplied treatment effect estimate
#' @param tau_hat_1 user-supplied treatment effect estimate for the setting with outcome added 1
#' @param mean_err error tolerance on the mean difference bewteen tau_hat and tau_hat_1
#'
#' @export
invariate_add_tests = function(tau_hat, tau_hat_1, mean_err=0.15) {
  print(abs(mean(tau_hat - tau_hat_1)))
  testthat::expect_equal(abs(mean(tau_hat - tau_hat_1)) < mean_err, TRUE)
}

#' @title helper function for testing treatment effect is invariant with a factor of 2 when outcome is multiplied with 2
#' @description helper function to test treatment effect is invariant with a factor of 2 when outcome is multiplied with 2
#' @param tau_hat user-supplied treatment effect estimate
#' @param tau_hat_2 user-supplied treatment effect estimate for the setting with outcome is multiplied with 2
#' @param mean_err error tolerance on the mean between 2x tau_hat and tau_hat_2
#'
#' @export
invariate_mult_tests = function(tau_hat, tau_hat_2, mean_err = 0.1) {
  print(abs(mean(2*tau_hat - tau_hat_2)) )
  testthat::expect_equal(abs(mean(2*tau_hat - tau_hat_2)) < mean_err, TRUE)
}
