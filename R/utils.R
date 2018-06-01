
#' For thresholding propensity scores
trim = function(x, min, max) {
	x[x>max] = max
	x[x<min] = min
	return(x)
}
