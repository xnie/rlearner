library(grf)
library(BART)

bart.cate = function(X, Y, W, X.test) {
    X.all = cbind(W, X)
    X.dummy = cbind(c(rep(0, nrow(X.test)), rep(1, nrow(X.test))), rbind(X.test, X.test))
    bart.fit = wbart(X.all, Y, X.dummy)
    bart.preds = colMeans(bart.fit$yhat.test)
    bart.tau = bart.preds[nrow(X.test) + 1:nrow(X.test)] - bart.preds[1:nrow(X.test)]
    bart.tau
}

bart.cate.crossfit = function(X, Y, W, X.test, nfold = 2) {
    X.all = cbind(W, X)
    test.preds = rep(0, nrow(X.test))
    train.preds = rep(NA, nrow(X))
    
    fold = sample.int(nfold, nrow(X), replace = TRUE)
    
    for (foldid in 1:nfold) {
        idx.out = which(fold == foldid)
        idx.in = which(fold != foldid)
        
        X.to.predict = rbind(X.test, X[idx.out,])
        X.dummy = cbind(c(rep(0, nrow(X.to.predict)),
                          rep(1, nrow(X.to.predict))),
                        rbind(X.to.predict, X.to.predict))
        bart.fit = wbart(X.all[idx.in,], Y[idx.in], X.dummy)
        bart.preds = colMeans(bart.fit$yhat.test)
        bart.tau.preds = bart.preds[nrow(X.to.predict) + 1:nrow(X.to.predict)] -
            bart.preds[1:nrow(X.to.predict)]
        
        test.preds = test.preds + bart.tau.preds[1:nrow(X.test)] / nfold
        train.preds[idx.out] = bart.tau.preds[nrow(X.test) + 1:length(idx.out)]
    }
    
    list(test.preds=test.preds, train.preds=train.preds)
}

nnls = function(M, v, constrained) {
    Dmat = t(M) %*% M
    dvec = t(M) %*% v
    Amat = matrix(0, ncol(M), sum(constrained))
    cons = which(constrained)
    for(iter in 1:length(cons)) {
        Amat[cons[iter], iter] = 1
    }
    bvec = rep(0, sum(constrained))
    soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec) 
    soln$solution
}