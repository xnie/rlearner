rm(list = ls())
library(ggplot2)

log = TRUE # log-log plots

tab = read.csv("output_boost.csv", header=TRUE)
algs = c("S", "T", "X", "U", "R", "RC", "causalboost", "oracle")
setups = c('A', 'B', 'C', 'D', 'E', 'F')
num.per.setup = 16

cols=c(1:7)

for (i in seq_along(setups)){
  data = tab[(1 + (i-1)*16):(i*16),algs]
  colnames(data) = c("s", "t", "x", "u", "r", "rc", "cb", "oracle")
  lograt = log(data[,1:7]/data[,8])

  pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  plot(1:2,1:2,type="n", xaxt="n", yaxt="n")
  legend <- legend("topleft", c("S", "T", "X", "U", "R", "RC", "CB"), pch = 7:13, col = cols, cex = 1.5, lwd = 2, lty = 0, plot = FALSE)

  if (log){
    pdf(paste("plots_boost/log-lasso", setups[i], ".pdf", sep=""))
  }
  else{
    pdf(paste("plots_boost/lasso", setups[i], ".pdf", sep=""))
  }


  if(log){
    ylim = range(c(log(unlist(data[,1:7]))))
    ylim[2] <- ylim[2]+legend$rect$h
    #if (i==4){
    #  ylim[2] <- ylim[2] + 2.6
    #}
    #if (setup==5){
    #  ylim[2] <- ylim[2] + 0.5
    #}
    plot(NA, NA, xlim = range(log(data[,8])),
         ylim=ylim,
         xlab = "oracle log-MSE",
         ylab = "log-MSE")
  }
  else{
    ylim = range(c(unlist(data[,1:7])))
    ylim[2] <- ylim[2]+legend$rect$h
    if (setup==4){
      ylim[2] <- ylim[2] + 2.6
    }
    if (setup==5){
      ylim[2] <- ylim[2] + 0.5
    }
    plot(NA, NA, xlim = range(data[,8]),
         ylim=ylim,
         xlab = "oracle MSE",
         ylab = "MSE")
  }
  for(iter in 1:length(algs)) {
    if (log){
      points(log(data[,7]), log(data[,iter]), col = cols[iter], pch = 6 + iter, lwd = 2, cex = 1.5)
    }
    else{
      points(data[,7], data[,iter], col = cols[iter], pch = 6 + iter, lwd = 2, cex = 1.5)
    }
  }
  legend <- legend("topleft", c("S", "T", "X", "U", "R", "RC", "CB"), pch = 7:13, col = cols, cex = 1.5, lwd = 2, lty = 0)
  abline(0, 1, lwd = 4, lty = 2)
  par = pardef
  dev.off()
}
