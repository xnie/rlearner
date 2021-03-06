% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skern.R
\name{skern}
\alias{skern}
\title{S-learner, implemented via kernel ridge regression with a Gaussian kernel}
\usage{
skern(
  x,
  w,
  y,
  k_folds = NULL,
  b_range = 10^(seq(-3, 3, 0.5)),
  lambda_range = 10^(seq(-3, 3, 0.5))
)
}
\arguments{
\item{x}{the input features}

\item{w}{the treatment variable (0 or 1)}

\item{y}{the observed response (real valued)}

\item{k_folds}{number of folds for cross validation}

\item{b_range}{the range of Gaussian kernel bandwidths for cross validation}

\item{lambda_range}{the range of ridge regression penalty factor for cross validation}
}
\value{
an skern object
}
\description{
S-learner, as proposed by Imai and Ratkovic (2013), implemented via
}
\examples{
\dontrun{
n = 100; p = 10

x = matrix(rnorm(n*p), n, p)
w = rbinom(n, 1, 0.5)
y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)

skern_fit = skern(x, w, y)
skern_est = predict(skern_fit, x)
}
}
