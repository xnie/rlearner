rm(list = ls())

library(rlearner)

args=(commandArgs(TRUE))
setup = as.numeric(args[1])
n = as.numeric(args[2])
p = as.numeric(args[3])
sigma = as.numeric(args[4])
alg = as.character(args[5])
NREP = as.numeric(args[6])
#
#setup = 2
#n=500
#p=6
#sigma=2
#alg='U'
#NREP=5


if (setup == 1) {

    get.params = function() {
        X = matrix(runif(n*p, min=0, max=1), n, p)
        b = sin(pi * X[,1] * X[,2]) + 2 * (X[,3] - 0.5)^2 + X[,4] + 0.5 * X[,5]
        e = sin(pi * X[,1] * X[,2])
        tau = (X[,1] + X[,2]) / 2
        list(X=X, b=b, tau=tau, e=e)
    }

} else if (setup == 2) {

    get.params = function() {
        X = matrix(rnorm(n*p), n, p)
        k=3
        rowm = rowMeans(X[,1:k] * sqrt(k))
        b = pmax(0, rowm)
        eta = 0.1
        e = pmax(eta, pmin(0.5 * (1 + sign(rowm) * rowm^2), 1-eta))
        tau = sin(2 * pi * X[,1])
        list(X=X, b=b, tau=tau, e=e)
    }

} else if (setup == 3) { # RCT

    get.params = function() {
        X = matrix(rnorm(n * p), n, p)
        b = pmax(0, X[,1] + X[,2], X[,3]) + pmax(0, X[,4] + X[,5])
        e = 0.5
        tau = X[,1] + log(1 + exp(X[,2]))
        list(X=X, b=b, tau=tau, e=e)
    }

} else if (setup == 4) { # constant treatment effect

    get.params = function() {
      X = matrix(rnorm(n * p), n, p)
      b = 2 * log(1 + exp(X[,1] + X[,2] + X[,3]))
      e = 1/(1 + exp(X[,2] + X[,3]))
      tau = rep(1, n)
      list(X=X, b=b, tau=tau, e=e)
    }

} else if (setup == 5) { # treat/control imbalance; complicated baseline+ treatment

    get.params = function() {
        X = matrix(rnorm(n * p), n, p)
        b = sin(pi * X[,1] * X[,2]) + (X[,3] + X[,4])^2
        e = 0.2
        tau = log(1 + exp(X[,3] + X[,5]))
        list(X=X, b=b, tau=tau, e=e)
    }

} else if (setup == 6) { # easy e, complicated baseline

    get.params = function() {
        X = matrix(rnorm(n*p), n, p)
        b = sqrt(pmax(0, X[,1] + X[,2], X[,3]))
        e = 1/(1 + exp(-X[,2]))
        tau = (X[,1] + X[,2])/2
        list(X=X, b=b, tau=tau, e=e)
    }

} else if (setup == 7) { # T

    get.params = function() {
        X = matrix(rnorm(n*p), n, p)
        b = (pmax(X[,1] + X[,2] + X[,3], 0) + pmax(X[,4] + X[,5], 0)) / 2
        e = 1/(1 + exp(-X[,1]) + exp(-X[,2]))
        tau = pmax(X[,1] + X[,2] + X[,3], 0) - pmax(X[,4] + X[,5], 0)
        list(X=X, b=b, tau=tau, e=e)
    }

} else if (setup == 8) { # setup 1 originally from the first paper draft

    get.params = function() {
        X = matrix(runif(n*p, min=0, max=1), n, p)
        b = 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
        e = sin(pi * X[,1] * X[,2])
        tau = (X[,1] + X[,2]) / 2
        list(X=X, b=b, tau=tau, e=e)
    }

} else {

    stop("bad setup")

}
    #params = get.params()
    #W = rbern(n, params$e)
    #Y = params$b + (W - 0.5) * params$tau + sigma * rnorm(n)
    #X.ns = do.call(cbind, lapply(1:p, function(col){matrix(ns(params$X[,col],df=7),n,7)}))

results.list = lapply(1:NREP, function(iter) {
    params = get.params()
    W = Rlab::rbern(n, params$e)
    Y = params$b + (W - 0.5) * params$tau + sigma * rnorm(n)
    X.ns = do.call(cbind, lapply(1:p, function(col){matrix(splines::ns(params$X[,col],df=7),n,7)}))

    if (alg == 'R') {

        est <- rlasso(X.ns, Y, W, lambda.choice = "lambda.min", constant.effect = TRUE, rs=FALSE)

    } else if (alg == 'RS') {

        est <- rlasso(X.ns, Y, W, lambda.choice = "lambda.min", constant.effect = TRUE, rs=TRUE)

    } else if (alg == 'Rnc') {

        est <- rlasso(X.ns, Y, W, lambda.choice = "lambda.min", constant.effect = FALSE, rs=FALSE)

    } else if (alg == 'RSnc') {

        est <- rlasso(X.ns, Y, W, lambda.choice = "lambda.min", constant.effect = FALSE, rs=TRUE)

    } else if (alg == 'S') {

        est <- slasso(X.ns, Y, W, lambda.choice = "lambda.min", constant.effect = TRUE)

    } else if (alg == 'T') {

        est <- tlasso(X.ns, Y, W, lambda.choice = "lambda.min")

    } else if (alg == 'X') {

        xlasso.fit <- xlasso(X.ns, Y, W, lambda.choice = "lambda.min")
        est <- list(tau.hat=predict(xlasso.fit)) # TODO

    } else if (alg == 'U') {

        est <- ulasso(X.ns, Y, W, lambda.choice = "lambda.min", cutoff=0.05)

    } else {

        stop("bad alg input")

    }

    est.mse = mean((est$tau.hat - params$tau)^2)
    print(est.mse)
    return(est.mse)
})
results = unlist(results.list, use.names=FALSE)
print(mean(results))
print(sd(results))

fnm = paste("results/output", alg, setup, n, p, sigma, NREP, "full.csv", sep="-")
write.csv(results, file=fnm)
