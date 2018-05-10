rm(list = ls())

library(rlearner)

# args=(commandArgs(TRUE))
# setup = as.numeric(args[1])
# n = as.numeric(args[2])
# p = as.numeric(args[3])
# sigma = as.numeric(args[4])
# k = as.numeric(args[5])
# eta = as.numeric(args[6])
# alg = as.numeric(args[7])
# NREP = as.numeric(args[8])

setup = 1
n=500
p=6
sigma=0.1
k=2
eta=0.1
alg='R'
NREP=50
rs=TRUE


if (setup == 1) {

    get.params = function() {
        X = matrix(runif(n*p, min=0, max=1), n, p)
        b = 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
        e = sin(pi * X[,1] * X[,2])
        tau = (X[,1] + X[,2]) / 2
        list(X=X, b=b, tau=tau, e=e)
    }

} else if (setup == 2) {

    get.params = function() {
        X = matrix(rnorm(n*p), n, p)
        rowm = rowMeans(X[,1:k] * sqrt(k))
        b = pmax(0, rowm)
        e = pmax(eta, pmin(0.5 * (1 + sign(rowm) * rowm^2), 1-eta))
        tau = sin(2 * pi * X[,1])
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
    W = rbern(n, params$e)
    Y = params$b + (W - 0.5) * params$tau + sigma * rnorm(n)
    X.ns = do.call(cbind, lapply(1:p, function(col){matrix(ns(params$X[,col],df=7),n,7)}))

    if (alg == 'R') {

        est <- rlasso(X.ns, Y, W, lambda.choice = "lambda.min", constant.effect = FALSE, standardize = FALSE, rs=rs)

    } else if (alg == 'S') {

        est <- slasso(X.ns, Y, W, lambda.choice = "lambda.min", constant.effect = TRUE, standardize = FALSE)

    } else if (alg == 'T') {

        est <- tlasso(X.ns, Y, W, lambda.choice = "lambda.min")

    } else if (alg == 'X') {

        est <- xlasso(X.ns, Y, W, lambda.choice = "lambda.min")

    } else if (alg == 'U') {

        est <- ulasso(X.ns, Y, W, lambda.choice = "lambda.min", cutoff=0)

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

#fnm = paste("results/output", alg, setup, n, p, sigma, k, eta, NREP, "full.csv", sep="-")
#write.csv(results, file=fnm)
