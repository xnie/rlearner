rm(list = ls())

library(xtable)
library(data.table)
args=(commandArgs(TRUE))
learner = as.character(args[1])
learner = "lasso"

if (learner != "boost" & learner != "lasso") {
  stop("learner needs to be boost or lasso")
}

filenames = list.files("results", pattern="*", full.names=TRUE)

param.names = c("alg", "learner", "setup", "n", "p", "sigma")
setup.values = c('A', 'B', 'C', 'D')

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

raw<-raw[raw$learner==learner & (raw$setup %in% setup.values), ] # only look at boost or lasso results

options(stringsAsFactors = FALSE)

raw = dcast(setDT(raw), setup + n + p + sigma ~ alg, value.var=c("mean"))

raw = raw[order(as.numeric(raw$sigma)),]
raw = raw[order(as.numeric(raw$p)),]
raw = raw[order(as.numeric(raw$n)),]
raw = raw[order(as.character(raw$setup)),]
rownames(raw) = 1:nrow(raw)

if (learner == "boost") {
  algs = c("S", "T", "X", "U", "causalboost", "R", "oracle")
  algs.tex = c("S", "T", "X", "U", "CB", "R", "oracle")
} else if (learner == "lasso") {
  algs = c("S", "T", "X", "U", "R", "RS", "oracle")
  algs.tex = algs
} else {
  stop("learner needs to be boost or lasso")
}
cols = c(c("setup", "n", "p", "sigma"), algs)
raw <- raw[, ..cols]
colnames(raw)[colnames(raw)=="causalboost"] <- "CB"

raw.round = raw
for (col in (5:11)){
  raw.round[,col] <- round(as.numeric(unlist(raw[,..col])),2)
}
raw = data.frame(apply(raw, 1:2, as.character))
raw.round = data.frame(apply(raw.round, 1:2, as.character))

# write raw csv output file
write.csv(raw, file=paste0("output_", learner, ".csv"))

# get a dataframe for each setup
raw.by.setup = lapply(c(setup.values), function(x) raw.round[raw.round$setup==x, ])

# write to latex tables
for (i in seq_along(setup.values)){
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
  print(setup.values[i])
  print(tab.setup)
  if (learner == "boost") {
    xtab.setup = xtable(tab.setup, caption = paste0("\\tt Mean squared error (MSE) running \\texttt{xgboost} from Setup ", setup.values[i], ". Results are averaged across 200 runs and rounded to two decimal places."), align="ccccccccccc", label=paste0("table:setup",i))
  } else if  (learner == "lasso"){
    xtab.setup = xtable(tab.setup, caption = paste0("\\tt Mean squared error (MSE) running \\texttt{lasso} from Setup ", setup.values[i], ". Results are averaged across 500 runs and rounded to two decimal places."), align="ccccccccccc", label=paste0("table:setup",i))
  }
  names(xtab.setup) <- c(c('n','p','$\\sigma$'), algs.tex)
  print(xtab.setup, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, file = paste("tables", "/simulation_results_setup_", setup.values[i], "_", learner, ".tex", sep=""))
}
