#' @import magrittr
#' @import caret
#' @import zeallot
#' @import dplyr
#' @import tidyr
#' @import purrr

# For thresholding propensity scores
trim = function(x, min, max) {
	x[x>max] = max
	x[x<min] = min
	return(x)
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
	n = 500
	x = stats::model.matrix(~.-1, data.frame("covariate_1" = rnorm(n), "covariate_2"= rnorm(n))) 
	logit_p = (x %*% c(1,1))
	p = exp(logit_p)/(1+exp(logit_p))
	w_bool = rbinom(n,1,p)==1
	w = factor(as.factor(w_bool %>% ifelse("treated", "control")), c("treated", "control"))
	tau = (x %*% c(1,1))^2
	m = x %*% c(1,-3)
	mu1 = m + tau/2
	mu0 = m - tau/2
	y = (m + tau/2*(2*w_bool-1))[,1] + rnorm(n)
	list(x, w, y, p, m, mu0, mu1, tau)
}
