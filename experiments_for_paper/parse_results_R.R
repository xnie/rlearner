rm(list = ls())

library(xtable)
library(data.table)
filenames = list.files("results_nc_oracle", pattern="*", full.names=TRUE)
#algvals=('Rnc' 'RSnc' 'Roracle' 'RSoracle' 'Rncoracle' 'RSncoracle')
param.names = c("alg", "setup", "n", "p", "sigma")
setup.values = c(1,2,3,4,5,6,7,8)

raw = data.frame(t(sapply(filenames, function(fnm) {

  output = read.csv(fnm)[,-1]
  params = strsplit(fnm, "-")[[1]][2:6]

  mse.mean = mean(output)
  mse.sd = sd(output)

  c(params,
    mse=sprintf("%.2f", round(mse.mean, 2)))
    #sd=sprintf("%.2f", round(mse.sd, 2)))
})))


rownames(raw) = 1:nrow(raw)


names(raw) = c(param.names,
               "mean")

options(stringsAsFactors = FALSE)
raw = data.frame(apply(raw, 1:2, as.character))

#raw = raw[!duplicated(raw[,1:5])] # can remove this line later
raw = dcast(setDT(raw), setup + n + p + sigma ~ alg, value.var=c("mean"))

raw = raw[order(as.numeric(raw$sigma)),]
raw = raw[order(as.numeric(raw$p)),]
raw = raw[order(as.numeric(raw$n)),]
raw = raw[order(as.numeric(raw$setup)),]
rownames(raw) = 1:nrow(raw)

raw = raw[,c('setup', 'n','p', 'sigma','R','Rnc', 'RS','RSnc','Roracle','Rncoracle','RSoracle','RSncoracle')]
write.csv(raw, file="output_R.csv")

# get a dataframe for each setup
raw.by.setup = lapply(c(setup.values), function(x) raw[raw$setup==x, ])

for (i in setup.values){
  tab.setup = cbind("", raw.by.setup[[i]][,-1])
  mse.idx = 1 + c(4:7)
  mse.oracle.idx = 1 + c(8:11)
  for(iter in 1:nrow(tab.setup)) {
    best.idx = mse.idx[which(as.numeric(tab.setup[iter,..mse.idx]) == min(as.numeric(tab.setup[iter,..mse.idx])))]
    best.oracle.idx = mse.oracle.idx[which(as.numeric(tab.setup[iter,..mse.oracle.idx]) == min(as.numeric(tab.setup[iter,..mse.oracle.idx])))]
    for (j in 1:length(best.idx)){
      best.idx.j = best.idx[j]
      tab.setup[iter,best.idx.j] = paste("\\bf", tab.setup[iter,..best.idx.j])
    }
    for (j in 1:length(best.oracle.idx)){
      best.oracle.idx.j = best.oracle.idx[j]
      tab.setup[iter,best.oracle.idx.j] = paste("\\bf", tab.setup[iter,..best.oracle.idx.j])
    }
  }
  tab.setup = tab.setup[,-1]
  print(i)
  print(tab.setup)
  xtab.setup = xtable(tab.setup, caption = paste("\\tt setup ", i, sep=""))
  names(xtab.setup) <- c('n','p','sigma', 'R', 'Rnc', 'RS', 'RSnc', 'Ro', 'Rnco', 'RSo', 'RSnco')
  print(xtab.setup, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, file = paste("tablesR/simulation_results_setup_", i, ".tex", sep=""))
}
