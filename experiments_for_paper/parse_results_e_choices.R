rm(list = ls())

library(xtable)
library(data.table)
filenames = list.files("results", pattern="*", full.names=TRUE)

param.names = c("alg", "setup", "n", "p", "sigma")
setup.values = c(1,2,3,4,5,6,7,8)
nms = c("S", "T", "X", "R", "RS") # add RSnc, Rnc, oracle

raw = data.frame(t(sapply(filenames, function(fnm) {

  output = read.csv(fnm)[,-1]
  params = strsplit(fnm, "-")[[1]][2:6]
  lambda.choice = strsplit(fnm, "-")[[1]][8]

  mse.mean = mean(output)
  lc = strsplit(lambda.choice, "\\.")[[1]][2]
  alg = paste(params[1], lc, sep='.')

  c(alg,
    params[2:5],
    mse=sprintf("%.8f", round(mse.mean, 8))) # change back to 2!
})))


rownames(raw) = 1:nrow(raw)


names(raw) = c(param.names,
               "mean")

options(stringsAsFactors = FALSE)

#raw = raw[!duplicated(raw[,1:5])] # can remove this line later
raw = dcast(setDT(raw), setup + n + p + sigma ~ alg, value.var=c("mean"))

raw = raw[order(as.numeric(raw$sigma)),]
raw = raw[order(as.numeric(raw$p)),]
raw = raw[order(as.numeric(raw$n)),]
raw = raw[order(as.numeric(raw$setup)),]
rownames(raw) = 1:nrow(raw)


raw = data.frame(apply(raw, 1:2, as.character))

write.csv(raw, file="output.csv")
#keep <- c("setup", "n", "p", "sigma", "RS1a.1se","RS1d.1se", "RSma.1se", "RSmd.1se", "RS.1se", "RSP.1se", "RS1a.min", "RS1d.min", "RSma.min", "RSmd.min", "RS.min", "RSP.min")
keep <- c("setup", "n", "p", "sigma", "R1a.1se","R1d.1se", "Rma.1se", "Rmd.1se", "R.1se", "R1a.min", "R1d.min", "Rma.min", "Rmd.min", "R.min")
raw <- raw[, (names(raw) %in% keep)]

logged = lapply(raw[,5:length(keep)], function(x) log(as.numeric(x)))
logged = do.call(cbind, lapply(1:(length(keep)-4), function(col)logged[[col]]))
total.row <- nrow(logged)
logged = logged[1:(total.row-16),]
colMeans(logged)

# get a dataframe for each setup
raw.by.setup = lapply(c(setup.values), function(x) raw[raw$setup==x, ])

for (i in setup.values){
  tab.setup = cbind("", raw.by.setup[[i]][,-1])
  mse.idx = 1 + c(4:11)
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
  xtab.setup = xtable(tab.setup, caption = paste("\\tt setup ", i, sep=""))
  names(xtab.setup) <- c('n','p','sig', "RS1a.1se","RS1d.1se", "RSma.1se", "RSmd.1se", "RS1a.min", "RS1d.min", "RSma.min", "RSmd.min")
  print(xtab.setup, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, file = paste("tables_RS_e_choices/simulation_results_setup_", i, ".tex", sep=""))
}
