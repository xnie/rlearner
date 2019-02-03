# The first file stacking_test has results for the first (continuous) tau specification from
# the paper. The file stacking_test_B has results for the second (discontinuous) one.


rm(list = ls())
library(RColorBrewer)

sigmas = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4)
in_use = 1:8

load("stacking_test.RData")

results_a = Reduce(rbind, lapply(all.results, rowMeans))
sd_a = Reduce(rbind, lapply(all.results, function(xx) apply(xx, 1, sd))) / sqrt(NREP)

load("stacking_test_B.RData")

results_b = Reduce(rbind, lapply(all.results, rowMeans))
sd_b = Reduce(rbind, lapply(all.results, function(xx) apply(xx, 1, sd))) / sqrt(NREP)

to_plot_a = data.frame(sigmas[in_use], sqrt(results_a[in_use, 2:4]))
names(to_plot_a) = c("sigma", "CF", "BART", "STACK")

to_plot_b = data.frame(sigmas[in_use], sqrt(results_b[in_use, 2:4]))
names(to_plot_b) = c("sigma", "CF", "BART", "STACK")

raw_a = sqrt(mean(results_a[,1]))
raw_b = sqrt(mean(results_b[,1]))

rng = range(c(to_plot_a[,2:4], to_plot_b[,2:4], raw_a, raw_b))

cols = brewer.pal(3, "Set1")

pdf("rstack_sharp.pdf")
pardef = par(xpd = FALSE, mar = c(4.5, 5, 3, 3) + 0.5, cex.lab=1.4, cex.axis=1.4, cex.main=1.4, cex.sub=1.4)
plot(to_plot_a$sigma, to_plot_a$CF, xlab = "sigma", ylab = "RMSE",
     log = "xy", ylim = rng, col = cols[1], pch = 15, cex = 1.75)
lines(to_plot_a$sigma, to_plot_a$CF, col = cols[1])
points(to_plot_a$sigma, to_plot_a$BART, col = cols[2], pch = 16, cex = 1.75)
lines(to_plot_a$sigma, to_plot_a$BART, col = cols[2])
points(to_plot_a$sigma, to_plot_a$STACK, col = cols[3], pch = 17, cex = 1.75)
lines(to_plot_a$sigma, to_plot_a$STACK, col = cols[3])
abline(h = raw_a, lty = 3, lwd = 2)
legend("bottomright", c("causal forest", "BART", "R-stack"),
       lty = 1, col = cols, pch = 15:17, pt.cex = 1.75, cex = 1.4)
par = pardef
dev.off()

pdf("rstack_smooth.pdf")
pardef = par(xpd = FALSE, mar = c(4.5, 5, 3, 3) + 0.5, cex.lab=1.4, cex.axis=1.4, cex.main=1.4, cex.sub=1.4)
plot(to_plot_b$sigma, to_plot_b$CF, xlab = "sigma", ylab = "RMSE",
     log = "xy", ylim = rng, col = cols[1], pch = 15, cex = 1.75)
lines(to_plot_b$sigma, to_plot_b$CF, col = cols[1])
points(to_plot_b$sigma, to_plot_b$BART, col = cols[2], pch = 16, cex = 1.75)
lines(to_plot_b$sigma, to_plot_b$BART, col = cols[2])
points(to_plot_b$sigma, to_plot_b$STACK, col = cols[3], pch = 17, cex = 1.75)
lines(to_plot_b$sigma, to_plot_b$STACK, col = cols[3])
abline(h = raw_b, lty = 3, lwd = 2)
legend("bottomright", c("causal forest", "BART", "R-stack"),
       lty = 1, col = cols, pch = 15:17, pt.cex = 1.75, cex = 1.4)
par = pardef
dev.off()
