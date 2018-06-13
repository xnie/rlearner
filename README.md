# rlearner: R-learner for Quasi-Oracle Estimation of Heterogenerous Treatment Effects

This package implements the R-learner for estimating
heterogeneous treatment effects, as proposed by Nie and Wager (2017). We consider a
setup where we observe data `(X, W, Y)` generated according
to the following general nonparameteric model

```
X ~ P(X)
W ~ P(W|X) where W is in {0,1}
Y = m(X) + 2(W-1)*tau(X) + e
```

With `e` being some mean-zero noise.

The R-learner estimates the heterogeneous treatment effect `tau(X)`.

To install this package in R, run the following commands:

```R
library(devtools) 
install_github("xnie/rlearner")
```
Example usage:

```R
library(rlearner)
library(zeallot)

# draw a sample of n observations from a simulation
data = toy_data_simulation(n) 
# data$x is a numeric matrix of covariates (each column having a name)
# data$w is a factor vector where the first level indicates "treatment" and the second "control"
# data$y is a numeric vector of outcomes

# Specify the machine learning models and 
# hyperparameters to be cross-validated over.
# This code specifies the use of elastic net
model_specs = list(
	glmnet = list(
	    tune_grid = expand.grid(
	       alpha=c(0,0.5,1),
	       lambda=exp(seq(-5,2,0.2))),
	    extra_args = list())
)

r.fit = rlearner_cv(data$x, data$w, data$y, model_specs=model_specs)
tau.hat = predict(r.fit, data$x)

print(paste("MSE of tau estimate:", mean((data$tau - tau.hat)^2), 2))
```

Using different machine learning models is as simple as changing the model specification:

```
# using gradient boosted tees
model_specs = list(
	gbm = list(
	    tune_grid = expand.grid(
	        n.trees = seq(1,501,20), 
	        interaction.depth=3, 
	        shrinkage = 0.1, 
	        n.minobsinnode=3),
	    extra_args = list(
	        verbose=F, 
	        bag.fraction=1))
)

r.fit = rlearner_cv(data$x, data$w, data$y, model_specs=model_specs)
tau.hat = predict(r.fit, data$x)
```

See `?learner_cv` for details of specifying the machine learning models and hyperparameters.

The package also implements S-, T-, U-, and X-learners. These can be called in a similar fashion:

```
t.fit = tlearner_cv(data$x, data$w, data$y, model_specs=model_specs)
tau.hat = predict(t.fit, data$x)
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
<b>Quasi-Oracle Estimation of Heterogeneous Treatment Effects.</b>
2017.
[<a href="https://arxiv.org/abs/1712.04912.pdf">arxiv</a>]
