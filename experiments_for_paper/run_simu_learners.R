rm(list = ls())

library(rlearner)
library(causalLearning)
library(magrittr)

start.time <- Sys.time()

#args=(commandArgs(TRUE))
#alg = as.character(args[1])
#learner = as.character(args[2])
#setup = as.character(args[3])
#n = as.numeric(args[4])
#p = as.numeric(args[5])
#sigma = as.numeric(args[6])
#NREP = as.numeric(args[7])
#
setup='C'
n=500
p=6
sigma=1
alg='R'
NREP=10
learner='boost'
print(alg)

if (setup == 'A') {
  get.params = function() {
    X = matrix(runif(n*p, min=0, max=1), n, p)
    b = sin(pi * X[,1] * X[,2]) + 2 * (X[,3] - 0.5)^2 + X[,4] + 0.5 * X[,5]
    eta = 0.1
    e = pmax(eta, pmin(sin(pi * X[,1] * X[,2]), 1-eta))
    tau = (X[,1] + X[,2]) / 2
    list(X=X, b=b, tau=tau, e=e)
  }

} else if (setup == 'B') { # RCT

  get.params = function() {
    X = matrix(rnorm(n * p), n, p)
    b = pmax(0, X[,1] + X[,2], X[,3]) + pmax(0, X[,4] + X[,5])
    e = 0.5
    tau = X[,1] + log(1 + exp(X[,2]))
    list(X=X, b=b, tau=tau, e=e)
  }

} else if (setup == 'C') { # constant treatment effect

  get.params = function() {
    X = matrix(rnorm(n * p), n, p)
    b = 2 * log(1 + exp(X[,1] + X[,2] + X[,3]))
    e = 1/(1 + exp(X[,2] + X[,3]))
    tau = rep(1, n)
    list(X=X, b=b, tau=tau, e=e)
  }

} else if (setup == 'D') { # T

  get.params = function() {
    X = matrix(rnorm(n*p), n, p)
    b = (pmax(X[,1] + X[,2] + X[,3], 0) + pmax(X[,4] + X[,5], 0)) / 2
    e = 1/(1 + exp(-X[,1]) + exp(-X[,2]))
    tau = pmax(X[,1] + X[,2] + X[,3], 0) - pmax(X[,4] + X[,5], 0)
    list(X=X, b=b, tau=tau, e=e)
  }
} else if (setup == 'E') {

  get.params = function() {
    X = matrix(rnorm(n*p), n, p)
    k=3
    rowm = rowMeans(X[,1:k] * sqrt(k))
    b = pmax(0, rowm)
    eta = 0.1
    e = pmax(eta, pmin(0.5 * (1 + sign(rowm) * rowm^2), 1-eta))
    tau = sin(2 * X[,1])
    list(X=X, b=b, tau=tau, e=e)
  }
} else if (setup == 'F') { # treat/control imbalance; complicated baseline+ treatment

  get.params = function() {
    X = matrix(rnorm(n * p), n, p)
    b = sin(pi * X[,1] * X[,2]) + (X[,3] + X[,4])^2
    e = 0.2
    tau = log(1 + exp(X[,3] + X[,5]))
    list(X=X, b=b, tau=tau, e=e)
  }

} else {

  stop("bad setup")

}
make_matrix = function(x) stats::model.matrix(~.-1, x)


results.list = lapply(1:NREP, function(iter) {

  params.train = get.params()
  W.train = rbinom(n,1,params.train$e)==1
  W.train.factor = factor(as.factor(W.train %>% ifelse("treated", "control")), c("treated", "control"))

  Y.train = params.train$b + (W.train - 0.5) * params.train$tau + sigma * rnorm(n)

  params.test = get.params()
  W.test = rbinom(n,1,params.test$e)==1
  W.test.factor = factor(as.factor(W.test %>% ifelse("treated", "control")), c("treated", "control"))
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

  if (learner == "lasso"){
    model_specs = list(
      glmnet = list(
        tune_grid = expand.grid(
          lambda=exp(seq(-5,2,0.2)),
          alpha=c(1)),
        extra_args = list()))
  } else if (learner == "boost"){
    model_specs = list(
      gbm = list(
          tune_grid = expand.grid(
              n.trees = seq(1,1001,20),
              interaction.depth=c(1,3),
              shrinkage = c(.05, .1),
              n.minobsinnode=c(3, 5, 10)),
          extra_args = list(
              verbose=F,
              bag.fraction=0.75))
      #xgbTree = list(
      #    tune_grid <-  expand.grid(eta = 0.1,
      #                            colsample_bytree=c(0.5,0.7),
      #                            max_depth=c(3,5),
      #                            nrounds=100,
      #                            gamma=1,
      #                            min_child_weight=c(2,5,8)),
      #                  extra_args=list())
    )
  }


  if (alg == 'R') {
    fitt = R_learner_cv(
      X.train, W.train.factor, Y.train,
      model_specs, model_specs, model_specs,
      economy=T)
    w.hat.oracle = params.train$e
    y.hat.oracle = params.train$b + (params.train$e-0.5) * params.train$tau
    #qplot(w.hat.oracle, fitt$p_hat)
    #qplot(y.hat.oracle, fitt$m_hat)

  } else if (alg == 'Ro') { # test if new R learner and old R learner implementation produce the same results

    fitt <- rlasso(X.train, Y.train, as.numeric(W.train), lambda.choice = "lambda.min", rs=FALSE)

  } else if (alg == "RC") {

    fitt = R_learner_cv(
      X.train, W.train.factor, Y.train,
      model_specs, model_specs, model_specs,
      economy=T, rc=TRUE)

    #w.hat.oracle = params.train$e
    #y.hat.oracle = params.train$b + (params.train$e-0.5) * params.train$tau
    #browser()
    #qplot(w.hat.oracle, fitt$p_hat)
    #qplot(y.hat.oracle, fitt$m_hat)

  } else if (alg == 'RS') {
    if (learner != "lasso"){
      stop("RS only has lasso implementation.")
    }

    fitt <- rlasso(X.train, Y.train, as.numeric(W.train), lambda.choice = "lambda.min", rs=TRUE)

  } else if (alg == 'oracleo') { # test if new R learner and old R learner implementation produce the same results

    w.hat.oracle = params.train$e
    y.hat.oracle = params.train$b + (params.train$e-0.5) * params.train$tau

    fitt <- rlasso(X.train, Y.train, as.numeric(W.train), lambda.choice = "lambda.min", rs=FALSE, w.hat = w.hat.oracle, y.hat = y.hat.oracle)

  } else if (alg == 'oracle') {

    w.hat.oracle = params.train$e
    y.hat.oracle = params.train$b + (params.train$e-0.5) * params.train$tau
    fitt = R_learner_cv(
      X.train, W.train.factor, Y.train,
      model_specs, model_specs, model_specs,
      m_hat = y.hat.oracle, p_hat = w.hat.oracle,
      economy=T)

  } else if (alg == 'S') {

    fitt = S_learner_cv(
      X.train, W.train.factor, Y.train,
      model_specs)

  } else if (alg == 'T') {

    fitt = T_learner_cv(
      X.train, W.train.factor, Y.train,
      model_specs)

  } else if (alg == 'X') {

    fitt = X_learner_cv(
      X.train, W.train.factor, Y.train,
      model_specs, model_specs, model_specs)

  } else if (alg == 'U') {

    fitt = U_learner_cv(
      X.train, W.train.factor, Y.train,
      model_specs, model_specs, model_specs,
      economy=T, p_min=0.05, p_max=0.95)

  } else if (alg == 'causalboost') {

    if (learner != "boost") {
      stop("causalboost is only available for learner=boost")
    }
    p_hat = xval_xfitt(X.train, W.train.factor, model_specs,
                      k_folds_cf=5, k_folds=5, economy=T, select_by="best")
    stratum = stratify(p_hat, W.train)$stratum
    fitt = cv.causalBoosting(X.train, as.numeric(W.train), Y.train, propensity=T, stratum=stratum)

  } else {

    stop("bad alg input")

  }
  tau.hat = predict(fitt, X.test)
  qplot(params.test$tau, tau.hat)
  #qplot(params.test$tau, fitt$tau.const)
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
