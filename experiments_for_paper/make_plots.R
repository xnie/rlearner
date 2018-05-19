rm(list = ls())
library(ggplot2)

log = TRUE# log-log plots

tab = read.csv("output.csv", header=TRUE)
algs = c("SP",  "T", "Xmd", "U", "Rmd", "RSmd", "oracle")
choice.ext = c(".1se", ".min")
setups = c(1:7)
num.per.setup = 16
ext.best <- rep(NA,length(algs))

cols=c(1:6)

for (i in 1:length(algs)){
  alg = algs[i]
  ext.best.alg = NA
  mse.best.alg = Inf
  for (ext in choice.ext) {
    data = tab[1:(num.per.setup*length(setups)), paste0(alg, ext)]
    logged = log(as.numeric(data))
    mse = mean(logged)
    if (mse < mse.best.alg){
      mse.best.alg = mse
      ext.best.alg = ext
    }
  }
  ext.best[[i]] = paste0(alg, ext.best.alg)
}
ext.best[[5]] = "Rmd.min"
ext.best[[6]] = "RSmd.min"


for (setup in setups){
  data = tab[(1 + (setup-1)*16):(setup*16),ext.best]
  colnames(data) = c("s", "t", "x", "u", "r", "rs", "oracle")
  lograt = log(data[,1:6]/data[,7])

  pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  plot(1:2,1:2,type="n", xaxt="n", yaxt="n")
  legend <- legend("topleft", c("S", "T", "X", "U", "R", "RS"), pch = 7:12, col = cols, cex = 1.5, lwd = 2, lty = 0, plot = FALSE)

  if (log){
    pdf(paste("plots/log-lasso", setup, ".pdf", sep=""))
  }
  else{
    pdf(paste("plots/lasso", setup, ".pdf", sep=""))
  }


  if(log){
    #ylim.maxmin = max(unlist(lapply(c(1:6), function(i) min(log(data[,i])))))
    ylim = range(c(log(unlist(data[,1:6]))))
    #ylim[2] <- ylim.maxmin+legend$rect$h
    ylim[2] <- ylim[2]+legend$rect$h
    if (setup==4){
      ylim[2] <- ylim[2] + 2.6
    }
    if (setup==5){
      ylim[2] <- ylim[2] + 0.5
    }
    plot(NA, NA, xlim = range(log(data[,7])),
         ylim=ylim,
         xlab = "oracle log-MSE",
         ylab = "log-MSE")
  }
  else{
    ylim = range(c(unlist(data[,1:6])))
    ylim[2] <- ylim[2]+legend$rect$h
    if (setup==4){
      ylim[2] <- ylim[2] + 2.6
    }
    if (setup==5){
      ylim[2] <- ylim[2] + 0.5
    }
    plot(NA, NA, xlim = range(data[,7]),
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
  legend <- legend("topleft", c("S", "T", "X", "U", "R", "RS"), pch = 7:12, col = cols, cex = 1.5, lwd = 2, lty = 0)
  abline(0, 1, lwd = 4, lty = 2)
  par = pardef
  dev.off()
}


#fvec1 = factor(c(rep("S", nrow(lograt1)),
#                 rep("T", nrow(lograt1)),
#                 rep("X", nrow(lograt1)),
#                 rep("R", nrow(lograt1))),
#               levels = c("S", "T", "X", "R"))
#
#df1 = data.frame(ratio=unlist(lograt1),
#                 label=fvec1)
#
#ggplot(df1, aes(x = label, y = ratio)) + geom_boxplot()
#
