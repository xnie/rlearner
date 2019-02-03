rm(list = ls())

library(foreach)

taufn = function(xx) { (xx[1] > 0) / (1 + exp(-xx[2])) }
mufn = function(xx) { 3 / (1 + exp(-xx[2] + xx[3])) }

prob = 0.5
p = 10
n = 10000

sigmas = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4)
NREP = 50

cl <- parallel::makeCluster(length(sigmas))
doParallel::registerDoParallel(cl)

all.results <- foreach(sigma = sigmas) %dopar% {
    
    source("lib.R")
    
    results = replicate(NREP, {
        
        X = matrix(rnorm(n*p), n, p)
        W = rbinom(n, 1, prob)
        mu0 = apply(X, 1, mufn)
        tau.train = apply(X, 1, taufn)
        Y = mu0 + W * tau.train + sigma * rnorm(n)
        
        n.test = n + 1
        X.test = matrix(rnorm(n.test*p), n.test, p)
        tau.test = apply(X.test, 1, taufn)
        
        cf = causal_forest(X, Y, W, W.hat = prob, tune.parameters = FALSE)
        tau.cf.oob = predict(cf)$predictions
        tau.cf = predict(cf, X.test)$predictions
        
        bout = capture.output(tau.bart <- bart.cate.crossfit(X, Y, W, X.test))

        RESP = Y - cf$Y.hat
        R.mat = cbind(1, W - prob,
                      (W - prob) * tau.cf.oob,
                      (W - prob) * tau.bart$train.preds)
        
        stack = nnls(R.mat, RESP, constrained = c(FALSE, FALSE, TRUE, TRUE))
        
        print("coefs")
        print(stack)
        tau.stack = stack[2] +
            stack[3] * tau.cf +
            stack[4] * tau.bart$test.preds
        
        out = c(var(tau.test),
                mean((tau.cf - tau.test)^2),
                mean((tau.bart$test.preds - tau.test)^2),
                mean((tau.stack - tau.test)^2))
        print("mse")
        print(out)
        out
    })
    
    print(paste("RESULTS FOR SIGMA", sigma))
    print(rowMeans(results))
    
    results
}

parallel::stopCluster(cl)

save.image("stacking_test.RData")
