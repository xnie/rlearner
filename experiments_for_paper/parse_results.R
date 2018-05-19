rm(list = ls())

library(xtable)
library(data.table)
filenames = list.files("results", pattern="*", full.names=TRUE)

param.names = c("alg", "setup", "n", "p", "sigma")
setup.values = c(1:7)

raw = data.frame(t(sapply(filenames, function(fnm) {

  output = read.csv(fnm)[,-1]
  params = strsplit(fnm, "-")[[1]][2:6]

  mse.mean = mean(output)

  c(params,
    mse=sprintf("%.8f", round(mse.mean, 8)))
})))

rownames(raw) = 1:nrow(raw)
names(raw) = c(param.names,
               "mean")

options(stringsAsFactors = FALSE)

raw = dcast(setDT(raw), setup + n + p + sigma ~ alg, value.var=c("mean"))

raw = raw[order(as.numeric(raw$sigma)),]
raw = raw[order(as.numeric(raw$p)),]
raw = raw[order(as.numeric(raw$n)),]
raw = raw[order(as.numeric(raw$setup)),]
rownames(raw) = 1:nrow(raw)
raw <- raw[,c("setup", "n", "p", "sigma", "S", "T", "X", "U", "R", "RS", "oracle")]

raw.round = raw
for (col in (5:11)){
  raw.round[,col] <- round(as.numeric(unlist(raw[,..col])),2)
}
raw = data.frame(apply(raw, 1:2, as.character))
raw.round = data.frame(apply(raw.round, 1:2, as.character))

# write raw csv output file
write.csv(raw, file="output.csv")

# get a dataframe for each setup
raw.by.setup = lapply(c(setup.values), function(x) raw.round[raw.round$setup==x, ])

# write to latex tables
for (i in setup.values){
  tab.setup = cbind("", raw.by.setup[[i]][,-1])
  mse.idx = 1 + c(4:9)
  for(iter in 1:nrow(tab.setup)) {
    best.idx = mse.idx[which(as.numeric(tab.setup[iter,mse.idx]) == min(as.numeric(tab.setup[iter,mse.idx])))]
    for (j in 1:length(best.idx)) {
      best.idx.j = best.idx[j]
      tab.setup[iter,best.idx.j] = paste("\\bf", tab.setup[iter,best.idx.j])
    }
  }
  tab.setup = tab.setup[,-1]
  print(i)
  print(tab.setup)
  xtab.setup = xtable(tab.setup, caption = paste0("\\tt Mean squared error (MSE) from Setup ", i, ". Results are averaged across 500 runs and rounded to two decimal places."), align="ccccccccccc", label=paste0("table:setup",i))
  names(xtab.setup) <- c('n','p','$\\sigma$', "S", "T", "X", "U", "R", "RS", "oracle")
  print(xtab.setup, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, file = paste("tables/simulation_results_setup_", i, ".tex", sep=""))
}
