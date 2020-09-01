#' @import magrittr
#' @import caret
#' @import zeallot
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import stringr

# For thresholding propensity scores
trim = function(x, min, max) {
	x[x>max] = max
	x[x<min] = min
	return(x)
}

sanitize_x = function(x){
	# make sure x is a numeric matrix with named columns (for caret)
	if (!is.matrix(x) | !is.numeric(x) | any(is.na(x))) {
		rlang::abort("x must be a numeric matrix with no missing values")
	}
	colnames(x) = str_c("covariate_", 1:ncol(x))
	return(x)
}

sanitize_input = function(x,w,y) {
  x = sanitize_x(x)

  if (!is.numeric(w)) {
		rlang::abort("the input w should be a numeric vector")
  }
	if (is.numeric(w) & all(w %in% c(0,1))) {
		w = w==1
	}

	# make sure y is a numeric vector
	if (!is.numeric(y)) {
		rlang::abort("y should be a numeric vector")
	}

	# make sure the dimensions align
	if (length(y)!=nrow(x) | length(w)!=nrow(x)) {
		rlang::abort("nrow(x), length(w), and length(y) should all be equal")
	}

	return(list(x,w,y))
}

#' @title Toy data simulation
#'
#' @description Generates a toy dataset of size \eqn{n} that can be used to experiment with the learners and meta-learners
#' in this package.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' toy_data_simulation # show the code- you can modify it to make your own simulations
#' library(zeallot) # imports the %<-% operator, which is syntactic sugar that performs multiple assignment out of a list
#' c(x, w, y, ...) %<-% toy_data_simulation(500) # draw a sample
#' # see what kind of objects these are
#' str(x)
#' str(w)
#' str(y)
#' @export
toy_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = rnorm(n), "covariate_2"= rnorm(n)))
	logit_p = (x %*% c(1,1))
	p = exp(logit_p)/(1+exp(logit_p))
	w = as.numeric(rbinom(n,1,p)==1)
	tau = (x %*% c(1,1))^2
	m = x %*% c(1,-3)
	mu1 = m + tau/2
	mu0 = m - tau/2
	y = (m + tau/2*(2*w-1))[,1] + rnorm(n)
	list(x=x, w=w, y=y, p=p, m=m, mu0=mu0, mu1=mu1, tau=tau)
}

#' @title Toy data simulation (easy mode)
#'
#' @description Generates a toy dataset of size \eqn{n} that can be used to experiment with the learners and meta-learners
#' in this package. The generative process should be easy to learn with linear methods.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' toy_data_simulation # show the code- you can modify it to make your own simulations
#' library(zeallot) # imports the %<-% operator, which is syntactic sugar that performs multiple assignment out of a list
#' c(x, w, y, ...) %<-% toy_data_simulation(500) # draw a sample
#' # see what kind of objects these are
#' str(x)
#' str(w)
#' str(y)
#' @export
easy_toy_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = rnorm(n), "covariate_2"= rnorm(n)))
	p = rep(0.5, n)
	w = as.numeric(rbinom(n,1,p)==1)
	tau = x %*% c(1,1)
	m = x %*% c(0.5,-0.5)
	mu1 = m + tau/2
	mu0 = m - tau/2
	y = (m + tau/2*(2*w-1))[,1]
	list(x=x, w=w, y=y, p=p, m=m, mu0=mu0, mu1=mu1, tau=tau)
}
#' @title Toy data simulation (easy mode) with continuous treatments
#'
#' @description Generates a toy dataset of size \eqn{n} that can be used to experiment with the learners and meta-learners
#' in this package. The generative process should be easy to learn with linear methods.
#'
#' @param n the number of samples to draw from the distribution
#' @return a list containing the covariate matrix, treatment vector, outcome vector, true propensity vector, true marginal outcome vector,
#' true control potential outcome vector, true treated potential outcome vector, and true treatment effect vector, in that order.
#' @examples
#' toy_data_simulation # show the code- you can modify it to make your own simulations
#' library(zeallot) # imports the %<-% operator, which is syntactic sugar that performs multiple assignment out of a list
#' c(x, w, y, ...) %<-% toy_data_simulation(500) # draw a sample
#' # see what kind of objects these are
#' str(x)
#' str(w)
#' str(y)
#' @export
continuous_toy_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = rnorm(n), "covariate_2"= rnorm(n)))
	w = runif(n, 0, 1)
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
#' toy_data_simulation # show the code- you can modify it to make your own simulations
#' library(zeallot) # imports the %<-% operator, which is syntactic sugar that performs multiple assignment out of a list
#' c(x, w, y, ...) %<-% t_toy_data_simulation(500) # draw a sample
#' # see what kind of objects these are
#' str(x)
#' str(w)
#' str(y)
#' @export
t_toy_data_simulation = function(n) {
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = rnorm(n), "covariate_2"= rnorm(n)))
	p = rep(0.5, n)
	w = as.numeric(rbinom(n,1,p)==1)
	mu1 = sin(x[,1] * 2)
	mu0 = x[,2] * 3 + 10
	y = w * mu1 + (1-w) * mu0
	tau = mu1 - mu0
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
  expect_equal(length(tau_hat), length(sim_data$tau))
  expect_equal(is.numeric(tau_hat), TRUE)
  learner_mse = mean((tau_hat - sim_data$tau)^2)
  print(learner_mse)
  expect_equal(learner_mse<mse, TRUE)
}
