
#' For thresholding propensity scores
threshold = function(x, min, max) {
	x[x>max] = max
	x[x<min] = min
	return(x)
}