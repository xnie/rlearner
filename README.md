# rlearner: R-learner for Quasi-Oracle Estimation of Hetereogenerous Treatment Effects

This package implements the R-learner for estimating
hetereogeneous treatment effects, as proposed by Nie and Wager (2017). We consider a
setup where we observe data `(X, Y, W)` generated according
to the following nonparameteric model,
```
Y = b(X) + (W-0.5) tau(X) + noise
```
and want to estimate the hetereogeneous treatment effect `tau(X)`.

To install this package in R, run the following commands:
```R
library(devtools) 
install_github("xnie/rlearner")
```
Example usage:

```R
library(rlearner)

n = 500; p = 12; sigma=1.0

X = matrix(rnorm(n * p), n, p)
b = sin(pi * X[,1] * X[,2]) + (X[,3] + X[,4])^2
e = 1/(1 + exp(-X[,2]))
tau = log(1 + exp(X[,3] + X[,5]))
W = rbinom(n, 1, e)
Y = b + (W - 0.5) * tau + sigma * rnorm(n)

r.fit = rlasso(X, Y, W, rs=FALSE)
tau.hat = predict(r.fit, newx = X)

print(paste("MSE of tau estimate:", round(mean((tau - tau.hat)^2), 2)))

```
To reproduce all simulation results from the paper (Nie and Wager, 2017), run the following commands:

```bash
cd experiments_for_paper
mkdir results
mkdir logging
mkdir tables
mkdir plots
./start.sh
source("parse_results.R")
source("make_plots.R")

```

#### References
Xinkun Nie and Stefan Wager.
<b>Quasi-Oracle Estimation of Hetereogeneous Treatment Effects.</b>
2017.
[<a href="https://arxiv.org/abs/1712.04912.pdf">arxiv</a>]
