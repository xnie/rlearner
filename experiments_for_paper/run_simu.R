rm(list = ls())

library(rlearner)
library(magrittr)

start_time <- Sys.time()

args=(commandArgs(TRUE))
alg = as.character(args[1])
learner = as.character(args[2])
setup = as.character(args[3])
n = as.numeric(args[4])
p = as.numeric(args[5])
sigma = as.numeric(args[6])
NREP = as.numeric(args[7])

#alg = "R"
#learner = "kernel"
#setup = "A"
#n = 500
#p = 6
#sigma = 2
#NREP = 20

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
results_list = lapply(1:NREP, function(iter) {
    print(paste0("starting iteration ", iter))

    params_train = get.params()
    W_train = rbinom(n, 1, params_train$e)
    Y_train = params_train$b + (W_train - 0.5) * params_train$tau + sigma * rnorm(n)

    params_test = get.params()
    W_test = rbinom(n, 1, params_test$e)
    Y_test = params_test$b + (W_test - 0.5) * params_test$tau + sigma * rnorm(n)

    while (T) {
      possibleError = tryCatch({ # in very rare occasions, there are numerical instability issues with the kernel ridge regression software (KLRS). In that case, we rerun it with the same simuated data.
    if (learner == "lasso") {
      X_ns = do.call(cbind, lapply(1:p, function(col){matrix(splines::ns(rbind(params_train$X, params_test$X)[,col],df=7), 2*n, 7)}))
      dim_ns = dim(X_ns)[2]
      X_ns = stats::model.matrix(~.*.-1, data.frame(X_ns)) # pairwise interaction (not including squared term for each column)
      X_ns_sq = do.call(cbind, lapply(1:dim_ns, function(col){matrix(X_ns[,col]^2)})) # squared term for each column
      X_ns = cbind(X_ns, X_ns_sq)
      X_train = data.frame(X_ns[1:n,]) %>% make_matrix
      X_test = data.frame(X_ns[(n+1):(2*n),]) %>% make_matrix
    } else if (learner == "boost" | learner == "kernel") {
      X_train = data.frame(params_train$X) %>% make_matrix
      X_test = data.frame(params_test$X) %>% make_matrix
    }
    else {
      stop("learner needs to be lasso, boost, or kernel.")
    }

    if (learner == "lasso") {
      if (alg == 'R') {

        fit <- rlasso(X_train, W_train, Y_train, lambda_choice = "lambda.min", rs = FALSE)

      } else if (alg == 'RS') {

        fit <- rlasso(X_train, W_train, Y_train, lambda_choice = "lambda.min", rs = TRUE)

      } else if (alg == 'oracle') {

        p_hat_oracle = params_train$e
        m_hat_oracle = params_train$b + (params_train$e-0.5) * params_train$tau
        fit <- rlasso(X_train, W_train, Y_train, lambda_choice = "lambda.min", rs = FALSE, p_hat = p_hat_oracle, m_hat = m_hat_oracle)

      } else if (alg == 'S') {

        fit <- slasso(X_train, W_train, Y_train, lambda_choice = "lambda.min", penalty_search=FALSE)

      } else if (alg == 'T') {

        fit <- tlasso(X_train, W_train, Y_train, lambda_choice = "lambda.min")

      } else if (alg == 'X') {

        fit <- xlasso(X_train, W_train, Y_train, lambda_choice = "lambda.min")

      } else if (alg == 'U') {

        fit <- ulasso(X_train, W_train, Y_train, lambda_choice = "lambda.1se", cutoff = 0.05)

      } else {

        stop("bad alg input")

      }
      tau_hat = predict(fit, newx=X_test)
    } else if (learner == "boost") {
      if (alg == 'R') {

        fit <- rboost(X_train, W_train, Y_train, nthread=1, verbose=TRUE)

      } else if (alg == 'oracle') {

        p_hat_oracle = params_train$e
        m_hat_oracle = params_train$b + (params_train$e-0.5) * params_train$tau
        fit <- rboost(X_train, W_train, Y_train, p_hat = p_hat_oracle, m_hat = m_hat_oracle, nthread = 1, verbose=TRUE)

      } else if (alg == 'S') {

        fit <- sboost(X_train, W_train, Y_train, nthread = 1, verbose = TRUE)

      } else if (alg == 'T') {

        fit <- tboost(X_train, W_train, Y_train, nthread = 1, verbose = TRUE)

      } else if (alg == 'X') {

        fit <- xboost(X_train, W_train, Y_train, nthread = 1, verbose = TRUE)

      } else if (alg == 'U') {

        fit <- uboost(X_train, W_train, Y_train, cutoff = 0.05, nthread = 1, verbose = TRUE)

      } else if (alg == 'causalboost') {
        causallearning_available = requireNamespace("causalLearning", quietly = TRUE)
        if (!causallearning_available) {
          stop("causalLearning needs to be installed for causal boosting.")
        }
        w_fit = cvboost(X_train, W_train, objective = "binary:logistic", nthread = 1, verbose = TRUE)
        p_hat = predict(w_fit)

        stratum = stratify(p_hat, W_train)$stratum
        fit = causalLearning::cv.causalBoosting(X_train, as.numeric(W_train), Y_train, propensity = T, stratum = stratum)

      } else {

        stop("bad alg input")

      }
      tau_hat = predict(fit, newx = X_test)
    }  else if (learner == "kernel"){
      krls_available = requireNamespace("KRLS2", quietly = TRUE)
      if (!krls_available) {
        stop("KRLS2 needs to be installed for kernel ridge regression.")
      } else{
        if (!packageVersion("KRLS2")>"1.1.0") {
        stop("KRLS2 needs to be of version at least 1.1.1.")
        }
      }

      if (alg == 'R') {

        fit = rkern(X_train, W_train, Y_train,
                         k_folds = 5,
                         b_range = 10^(seq(-3,3,0.5)),
                         lambda_range = 10^(seq(-3,3,0.5)))

      } else if (alg == 'oracle') {

        p_hat_oracle = params_train$e
        m_hat_oracle = params_train$b + (params_train$e-0.5) * params_train$tau

        fit = rkern(X_train, W_train, Y_train,
                    p_hat = p_hat_oracle,
                    m_hat = m_hat_oracle,
                    k_folds = 5,
                    b_range = 10^(seq(-3,3,0.5)),
                    lambda_range = 10^(seq(-3,3,0.5)))

      } else if (alg == 'S') {

        fit <- skern(X_train, W_train, Y_train,
                     k_folds = 5,
                     b_range = 10^(seq(-3,3,0.5)),
                     lambda_range = 10^(seq(-3,3,0.5)))
      }

      else if (alg == 'T') {

        fit <- tkern(X_train, W_train, Y_train,
                     k_folds = 5,
                     b_range = 10^(seq(-3,3,0.5)),
                     lambda_range = 10^(seq(-3,3,0.5)))

      }

      else if (alg == 'X') {

        fit <- xkern(X_train, W_train, Y_train,
                    k_folds=NULL,
                    b_range =10^(seq(-3,3,0.5)),
                    lambda_range = 10^(seq(-3,3,0.5)))

      }
      else if (alg == 'U') {

        fit <- ukern(X_train, W_train, Y_train,
                     k_folds = 5,
                     b_range = 10^(seq(-3,3,0.5)),
                     lambda_range = 10^(seq(-3,3,0.5)))
      }
      else {

        stop("bad alg input")

      }
      tau_hat = predict(fit, newx = X_test)

    }
    else {

      stop("the learner needs to be lasso, boost, or kernel")

    }

    est.mse = mean((tau_hat - params_test$tau)^2)
    end_time <- Sys.time()
    time_taken <- end_time - start_time
    print(time_taken)
    print(est.mse)
    },
    error = function(e) {e})
    if (!inherits(possibleError, "error")){
      if (est.mse < 1000) {
        break
      }
      else {
        print(paste0("moving on. mse is ", est.mse))
        next()
      }
    } else {
      print("moving on. error occured")
      next()
    }
  }
  return(est.mse)
})
results = unlist(results_list, use.names=FALSE)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

fnm = paste("results/output", alg, learner, setup, n, p, sigma, NREP, "full.csv", sep="-")
write.csv(results, file=fnm)
