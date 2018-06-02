rm(list = ls())

library(rlearner)

start.time <- Sys.time()

#args=(commandArgs(TRUE))
#alg = as.character(args[1])
#setup = as.numeric(args[2])
#n = as.numeric(args[3])
#p = as.numeric(args[4])
#sigma = as.numeric(args[5])
#NREP = as.numeric(args[6])
#
setup=1
n=500
p=6
sigma=0.1
alg='T'
NREP=10
print(alg)

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

} else {

    stop("bad setup")

}

results.list = lapply(1:NREP, function(iter) {

    params.train = get.params()
    W.train = Rlab::rbern(n, params.train$e)
    Y.train = params.train$b + (W.train - 0.5) * params.train$tau + sigma * rnorm(n)

    params.test = get.params()
    W.test = Rlab::rbern(n, params.test$e)
    Y.test = params.test$b + (W.test - 0.5) * params.test$tau + sigma * rnorm(n)

    X.ns = do.call(cbind, lapply(1:p, function(col){matrix(splines::ns(rbind(params.train$X, params.test$X)[,col],df=7), 2*n, 7)}))
    X.ns.train = X.ns[1:n,]
    X.ns.test = X.ns[(n+1):(2*n),]

    if (alg == 'R') {

        r.fit <- rlasso(X.ns.train, Y.train, W.train, lambda.choice = "lambda.min", rs=FALSE)
        tau.hat <- predict(r.fit, newx=X.ns.test)

    } else if (alg == 'RS') {

        rs.fit <- rlasso(X.ns.train, Y.train, W.train, lambda.choice = "lambda.min", rs=TRUE)
        tau.hat <- predict(rs.fit, newx=X.ns.test)

    } else if (alg == 'oracle') {

        w.hat.oracle = params.train$e
        y.hat.oracle = params.train$b + (params.train$e-0.5) * params.train$tau
        oracle.fit <- rlasso(X.ns.train, Y.train, W.train, lambda.choice = "lambda.min", rs=FALSE, w.hat=w.hat.oracle, y.hat=y.hat.oracle)
        tau.hat <- predict(oracle.fit, newx=X.ns.test)

    } else if (alg == 'S') {

        s.fit <- slasso(X.ns.train, Y.train, W.train, lambda.choice = "lambda.min", penalty.search=TRUE)
        tau.hat <- predict(s.fit, newx=X.ns.test)

    } else if (alg == 'T') {

        t.fit <- tlasso(X.ns.train, Y.train, W.train, lambda.choice = "lambda.min")
        tau.hat <- predict(t.fit, newx=X.ns.test, s="lambda.min")

    } else if (alg == 'X') {

        x.fit <- xlasso(X.ns.train, Y.train, W.train, lambda.choice = "lambda.min")
        tau.hat <- predict(x.fit, newx=X.ns.test, s="lambda.min")

    } else if (alg == 'U') {

        u.fit <- ulasso(X.ns.train, Y.train, W.train, lambda.choice = "lambda.1se", cutoff = 0.05)
        tau.hat <- predict(u.fit, newx=X.ns.test, s="lambda.1se")

    } else {

        stop("bad alg input")

    }

    est.mse = mean((tau.hat - params.test$tau)^2)
    print(est.mse)
    return(est.mse)
})
results = unlist(results.list, use.names=FALSE)
print(mean(results))
print(sd(results))
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

fnm = paste("results/output", alg, setup, n, p, sigma, NREP, "full.csv", sep="-")
write.csv(results, file=fnm)
