rm(list = ls())

library(xtable)
library(data.table)
filenames = list.files("results", pattern="*", full.names=TRUE)

param.names = c("alg", "learner", "setup", "n", "p", "sigma")
setup.values = c('A', 'B', 'C', 'D', 'E', 'F')

raw = data.frame(t(sapply(filenames, function(fnm) {

  output = read.csv(fnm)[,-1]
  params = strsplit(fnm, "-")[[1]][2:7]

  mse.mean = mean(output)

  c(params,
    mse=sprintf("%.8f", round(mse.mean, 8)))
})))


rownames(raw) = 1:nrow(raw)
names(raw) = c(param.names,
               "mean")

raw<-raw[raw$learner=="boost",] # only look at boost results

options(stringsAsFactors = FALSE)

raw = dcast(setDT(raw), setup + n + p + sigma ~ alg, value.var=c("mean"))

raw = raw[order(as.numeric(raw$sigma)),]
raw = raw[order(as.numeric(raw$p)),]
raw = raw[order(as.numeric(raw$n)),]
raw = raw[order(as.character(raw$setup)),]
rownames(raw) = 1:nrow(raw)
raw <- raw[,c("setup", "n", "p", "sigma", "S", "T", "X", "U", "causalboost", "R", "RC", "oracle")]

raw.round = raw
for (col in (5:12)){
  raw.round[,col] <- round(as.numeric(unlist(raw[,..col])),2)
}
raw = data.frame(apply(raw, 1:2, as.character))
raw.round = data.frame(apply(raw.round, 1:2, as.character))

# write raw csv output file
write.csv(raw, file="output_boost.csv")

# get a dataframe for each setup
raw.by.setup = lapply(c(setup.values), function(x) raw.round[raw.round$setup==x, ])

# write to latex tables
for (i in seq_along(setup.values)){
  tab.setup = cbind("", raw.by.setup[[i]][,-1])
  mse.idx = 1 + c(4:10)
  for(iter in 1:nrow(tab.setup)) {
    best.idx = mse.idx[which(as.numeric(tab.setup[iter,mse.idx]) == min(as.numeric(tab.setup[iter,mse.idx])))]
    for (j in 1:length(best.idx)) {
      best.idx.j = best.idx[j]
      tab.setup[iter,best.idx.j] = paste("\\bf", tab.setup[iter,best.idx.j])
    }
  }
  tab.setup = tab.setup[,-1]
  print(setup.values[i])
  print(tab.setup)
  xtab.setup = xtable(tab.setup, caption = paste0("\\tt Mean squared error (MSE) running gbm from Setup ", setup.values[i], ". Results are averaged across 500 runs and rounded to two decimal places."), align="cccccccccccc", label=paste0("table:setup",i))
  names(xtab.setup) <- c('n','p','$\\sigma$', "S", "T", "X", "U", "cb", "R", "RC", "oracle")
  print(xtab.setup, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, file = paste("tables_boost/simulation_results_setup_", setup.values[i], ".tex", sep=""))
}
