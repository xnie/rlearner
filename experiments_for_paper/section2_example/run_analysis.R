rm(list = ls())

set.seed(1)

library(rlearner)
library(glmnet)
library(dbarts)
library(RColorBrewer)

proc.time()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data = read.csv("data_clean.csv")

proc.time()

idx.all = sample(c(sample(which(data$W == 0), sum(data$W) * 3/2), which(data$W == 1)))

n.train = 100000
n.test = 25000
n.holdout = length(idx.all) - n.train - n.test

#
# MAKE SYNTHETIC TAU
#

X.all = as.matrix(data[idx.all,1:11])
Y.obs.all = data[idx.all,12]
W.all = data[idx.all,13]

TAU.all = - X.all[,"vote00"] * 0.5 / (1 + 50 / X.all[,"age"])

FLIP = rbinom(length(idx.all), 1, abs(TAU.all))
Y.PO = t(sapply(1:length(idx.all), function(ii) {
  if (FLIP[ii] == 0) {
    return(rep(Y.obs.all[ii], 2))
  } else if (TAU.all[ii] > 0) {
    return(c(0, 1))
  } else {
    return(c(1, 0))
  }
}))

plot(smooth.spline(TAU.all, Y.PO[,2] - Y.PO[,1]))
abline(0, 1)

Y.all = sapply(1:length(idx.all), function(ii) {
  Y.PO[ii, W.all[ii] + 1]
})

proc.time()

# Train data
X = X.all[1:n.train,]
W = W.all[1:n.train]
Y = Y.all[1:n.train]

# Test data
X.test = X.all[n.train + 1:n.test,]
W.test = W.all[n.train + 1:n.test]
Y.test = Y.all[n.train + 1:n.test]
TAU.test = TAU.all[n.train + 1:n.test]

# Holdout data
X.holdout = X.all[n.train + n.test + 1:n.holdout,]
W.holdout = W.all[n.train + n.test + 1:n.holdout]
Y.holdout = Y.all[n.train + n.test + 1:n.holdout]


#
# fit propensity model
#

proc.time()
W.boost = cvboost(X, W, objective = "binary:logistic", nthread = 4)
W.hat.boost = predict(W.boost)
proc.time()

W.lasso = cv.glmnet(X, W, family = "binomial", keep = TRUE)
W.hat.lasso =
  W.lasso$fit.preval[,!is.na(colSums(W.lasso$fit.preval))][, W.lasso$lambda == W.lasso$lambda.min]

# boosting wins
print(round(c(-mean(W * log(W.hat.boost) + (1 - W) * log(1 - W.hat.boost)),
              -mean(W * log(W.hat.lasso) + (1 - W) * log(1 - W.hat.lasso))), 4))
W.hat = W.hat.boost


#
# fit marginal response model
#

proc.time()
Y.boost = cvboost(X, Y, objective = "binary:logistic", nthread = 4)
Y.hat.boost = predict(Y.boost)
proc.time()

Y.lasso = cv.glmnet(X, Y, keep = TRUE, family = "binomial")
Y.hat.lasso =
  Y.lasso$fit.preval[,!is.na(colSums(Y.lasso$fit.preval))][, Y.lasso$lambda == Y.lasso$lambda.min]

# boosting wins
print(round(c(-mean(Y * log(Y.hat.boost) + (1 - Y) * log(1 - Y.hat.boost)),
              -mean(Y * log(Y.hat.lasso) + (1 - Y) * log(1 - Y.hat.lasso))), 4))
Y.hat = Y.hat.boost

#
# fit R-learner given chosen nuisance components
#

proc.time()
tau.boost = rboost(X, Y, W, w.hat = W.hat, y.hat = Y.hat, nthread = 4)
tau.hat.boost = predict(tau.boost)
proc.time()

tau.lasso = rlasso(X, Y, W, w.hat = W.hat, y.hat = Y.hat)
tau.hat.lasso = predict(tau.lasso)

# who wins on CV error?
print(round(c(mean((Y - Y.hat - tau.hat.boost * (W - W.hat))^2),
              mean((Y - Y.hat - tau.hat.lasso * (W - W.hat))^2)), 4))

# lasso wins on holdout
tau.hat.boost.holdout = predict(tau.boost, X.holdout)
tau.hat.lasso.holdout = predict(tau.lasso, X.holdout)

Y.hat.holdout = predict(Y.boost, X.holdout)
W.hat.holdout = predict(W.boost, X.holdout)

print(round(c(mean((Y.holdout - Y.hat.holdout - tau.hat.boost.holdout * (W.holdout - W.hat.holdout))^2),
              mean((Y.holdout - Y.hat.holdout - tau.hat.lasso.holdout * (W.holdout - W.hat.holdout))^2)), 4))

# lasso in fact wins on oracle test set error
tau.hat.boost.test = predict(tau.boost, X.test)
tau.hat.lasso.test = predict(tau.lasso, X.test)

print(round(c(mean((TAU.test - tau.hat.boost.test)^2),
              mean((TAU.test - tau.hat.lasso.test)^2)), 4))


#
# try the S-lasso
#

proc.time()
tau.slasso = slasso(X, Y, W)
tau.hat.slasso.test = predict(tau.slasso, X.test)
round(mean((TAU.test - tau.hat.slasso.test)^2), 4)
proc.time()

#
# try PS-BART
#

XpsW = cbind(X, W.hat, W)

W.hat.test = predict(W.boost, X.test)
Xps01 = rbind(cbind(X.test, W.hat.test, 0),
              cbind(X.test, W.hat.test, 1))

proc.time()

bart.fit = dbarts::bart(XpsW, Y, x.test = Xps01)

bart.post = colMeans(1 / (1 + exp(-bart.fit$yhat.test)))
mu0.hat.bart.test = bart.post[1:nrow(X.test)]
mu1.hat.bart.test = bart.post[nrow(X.test) + 1:nrow(X.test)]
tau.hat.bart.test = mu1.hat.bart.test - mu0.hat.bart.test

round(mean((TAU.test - tau.hat.bart.test)^2), 4)

proc.time()

#
# make summaries
#

print(round(c(mean((TAU.test - tau.hat.boost.test)^2),
              mean((TAU.test - tau.hat.lasso.test)^2),
              mean((TAU.test - tau.hat.slasso.test)^2),
              mean((TAU.test - tau.hat.bart.test)^2),
              var(TAU.test)), 5))


cols = brewer.pal(3, "Set2")

pdf("tau_histogram.pdf")
pardef = par(mar = c(4.5, 5, 3, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(TAU.all[1:n.train], main = "", xlab = "CATE", col = cols[1])
dev.off()

TAU.bucket = as.numeric(cut(TAU.test,
                            breaks = c(quantile(TAU.test[TAU.test < 0], c(0, 0.25, 0.5, 0.75, 1)) - 0.001, 0.001))) - 1
TAU.mids = sapply(0:4, function(ii) median(TAU.test[TAU.bucket == ii]))
ATmat = rbind(TAU.mids, TAU.mids)
ATmat = ATmat + c(-0.005, 0.005)
ATvec = c(ATmat)

BPDF = data.frame(CATE=c(2 * TAU.bucket, 1 + 2 * TAU.bucket),
                  prediction=c(tau.hat.lasso.test, tau.hat.boost.test))

pdf("tau_boxplot.pdf")
pardef = par(xpd = FALSE, mar = c(4.5, 5, 3, 8) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
boxplot(prediction ~ CATE, data = BPDF, col = cols[2:3],
        at = ATvec,  pars = list(boxwex = 0.007), xlim = range(ATvec),
        names = rep("", 10), xlab = "CATE", ylab = "prediction", xaxt = "n")
axis(1, at = round(TAU.mids, 2))
abline(0, 1, lwd = 2)

par(xpd = TRUE)
legend(0.03, 0.1, c("lasso", "boost"), fill = cols[2:3], cex = 1.5)
par = pardef
dev.off()

save.image("analysis_results.RData")

proc.time()
