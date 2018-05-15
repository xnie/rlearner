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

#logged = lapply(raw[,5:18], function(x) log(as.numeric(x)))
#logged = do.call(cbind, lapply(1:14, function(col)logged[[col]]))
#total.row <- nrow(logged)
#logged = logged[1:(total.row-16),]
#colMeans(logged)

raw = data.frame(apply(raw, 1:2, as.character))

write.csv(raw, file="output.csv")

# get a dataframe for each setup
raw.by.setup = lapply(c(setup.values), function(x) raw[raw$setup==x, ])

for (i in setup.values){
  tab.setup = cbind("", raw.by.setup[[i]][,-1])
  mse.idx = 1 + c(4:9)
  for (alg.idx in 1:7){
    mse.idx=(5 + (alg.idx-1)*2):(6+(alg.idx-1)*2)
    for(iter in 1:nrow(tab.setup)) {
      best.idx = mse.idx[which(as.numeric(tab.setup[iter,..mse.idx]) == min(as.numeric(tab.setup[iter,..mse.idx])))]
      for (j in 1:length(best.idx)){
        best.idx.j = best.idx[j]
        tab.setup[iter,best.idx.j] = paste("\\bf", tab.setup[iter,..best.idx.j])
      }
    }
  }
  tab.setup = tab.setup[,-1]
  print(i)
  print(tab.setup)
  xtab.setup = xtable(tab.setup, caption = paste("\\tt setup ", i, sep=""))
  names(xtab.setup) <- c('n','p','sig', 'R.1', 'R.m', 'RS.1', 'RS.m', 'S.1', 'S.m', 'T.1', 'T.m', 'U.1', 'U.m', 'X.1', 'X.m', 'o.1', 'o.m' )
  print(xtab.setup, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, file = paste("tables_lambdachoice/simulation_results_setup_", i, ".tex", sep=""))
}
