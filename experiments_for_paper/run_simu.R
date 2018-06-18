rm(list = ls())

library(rlearner)
library(causalLearning)
library(xgboost)
library(magrittr)

start.time <- Sys.time()

args=(commandArgs(TRUE))
alg = as.character(args[1])
learner = as.character(args[2])
setup = as.character(args[3])
n = as.numeric(args[4])
p = as.numeric(args[5])
sigma = as.numeric(args[6])
NREP = as.numeric(args[7])

if (setup == 'A') {
  get.params = function() {
    X = matrix(runif(n*p, min=0, max=1), n, p)
    b = sin(pi * X[,1] * X[,2]) + 2 * (X[,3] - 0.5)^2 + X[,4] + 0.5 * X[,5]
    eta = 0.1
    e = pmax(eta, pmin(sin(pi * X[,1] * X[,2]), 1-eta))
    tau = (X[,1] + X[,2]) / 2
    list(X=X, b=b, tau=tau, e=e)
  }

} else if (setup == 'B') {

  get.params = function() {
    X = matrix(rnorm(n * p), n, p)
    b = pmax(0, X[,1] + X[,2], X[,3]) + pmax(0, X[,4] + X[,5])
    e = 0.5
    tau = X[,1] + log(1 + exp(X[,2]))
    list(X=X, b=b, tau=tau, e=e)
  }

} else if (setup == 'C') {

  get.params = function() {
    X = matrix(rnorm(n * p), n, p)
    b = 2 * log(1 + exp(X[,1] + X[,2] + X[,3]))
    e = 1/(1 + exp(X[,2] + X[,3]))
    tau = rep(1, n)
    list(X=X, b=b, tau=tau, e=e)
  }

} else if (setup == 'D') {

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

make_matrix = function(x) stats::model.matrix(~.-1, x)

results.list = lapply(1:NREP, function(iter) {

    params.train = get.params()
    W.train = rbinom(n, 1, params.train$e)
    Y.train = params.train$b + (W.train - 0.5) * params.train$tau + sigma * rnorm(n)

    params.test = get.params()
    W.test = rbinom(n, 1, params.test$e)
    Y.test = params.test$b + (W.test - 0.5) * params.test$tau + sigma * rnorm(n)

    if (learner == "lasso") {
      X.ns = do.call(cbind, lapply(1:p, function(col){matrix(splines::ns(rbind(params.train$X, params.test$X)[,col],df=7), 2*n, 7)}))
      dim.ns = dim(X.ns)[2]
      X.ns = stats::model.matrix(~.*.-1, data.frame(X.ns)) # pairwise interaction (not including squared term for each column)
      X.ns.sq = do.call(cbind, lapply(1:dim.ns, function(col){matrix(X.ns[,col]^2)})) # squared term for each column
      X.ns = cbind(X.ns, X.ns.sq)
      X.train = data.frame(X.ns[1:n,]) %>% make_matrix
      X.test = data.frame(X.ns[(n+1):(2*n),]) %>% make_matrix
    }
    else if (learner == "boost") {
      X.train = data.frame(params.train$X) %>% make_matrix
      X.test = data.frame(params.test$X) %>% make_matrix
    }
    else {
      stop("learner needs to be lasso or boost.")
    }

    if (learner == "lasso") {
      if (alg == 'R') {

        fit <- rlasso(X.train, Y.train, W.train, lambda.choice = "lambda.min", rs=FALSE)

      } else if (alg == 'RS') {

        fit <- rlasso(X.train, Y.train, W.train, lambda.choice = "lambda.min", rs=TRUE)

      } else if (alg == 'oracle') {

        w.hat.oracle = params.train$e
        y.hat.oracle = params.train$b + (params.train$e-0.5) * params.train$tau
        fit <- rlasso(X.train, Y.train, W.train, lambda.choice = "lambda.min", rs=FALSE, w.hat=w.hat.oracle, y.hat=y.hat.oracle)

      } else if (alg == 'S') {

        fit <- slasso(X.train, Y.train, W.train, lambda.choice = "lambda.min", penalty.search=TRUE)

      } else if (alg == 'T') {

        fit <- tlasso(X.train, Y.train, W.train, lambda.choice = "lambda.min")

      } else if (alg == 'X') {

        fit <- xlasso(X.train, Y.train, W.train, lambda.choice = "lambda.min")

      } else if (alg == 'U') {

        fit <- ulasso(X.train, Y.train, W.train, lambda.choice = "lambda.1se", cutoff = 0.05)

      } else {

        stop("bad alg input")

      }
      tau.hat = predict(fit, newx=X.test)
    }

    else if (learner == "boost") {
      if (alg == 'R') {

        fit <- rboost(X.train, Y.train, W.train, nthread=1)

      } else if (alg == 'oracle') {

        w.hat.oracle = params.train$e
        y.hat.oracle = params.train$b + (params.train$e-0.5) * params.train$tau
        fit <- rboost(X.train, Y.train, W.train, w.hat=w.hat.oracle, y.hat=y.hat.oracle, nthread=1)

      } else if (alg == 'S') {

        fit <- sboost(X.train, Y.train, W.train, nthread=1)

      } else if (alg == 'T') {

        fit <- tboost(X.train, Y.train, W.train, nthread=1)

      } else if (alg == 'X') {

        fit <- xboost(X.train, Y.train, W.train, nthread=1)

      } else if (alg == 'U') {

        fit <- uboost(X.train, Y.train, W.train, cutoff=0.05, nthread=1)

      } else if (alg == 'causalboost') {

        w.fit = cvboost(X.train, W.train, objective="binary:logistic", nthread=1)
        w.hat = predict(w.fit)

        stratum = stratify(w.hat, W.train)$stratum
        fit = cv.causalBoosting(X.train, as.numeric(W.train), Y.train, propensity=T, stratum=stratum)

      } else {

        stop("bad alg input")

      }
      tau.hat = predict(fit, newx=X.test)
    }
    else {

      stop("the learner needs to be lasso or boost")

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

fnm = paste("results/output", alg, learner, setup, n, p, sigma, NREP, "full.csv", sep="-")
write.csv(results, file=fnm)
