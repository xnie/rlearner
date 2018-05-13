rm(list = ls())

library(xtable)
library(data.table)
filenames = list.files("results", pattern="*", full.names=TRUE)

param.names = c("alg", "setup", "n", "p", "sigma")
nms = c("S", "T", "X", "R", "RS") # add RSnc, Rnc, oracle

raw = data.frame(t(sapply(filenames, function(fnm) {

	output = read.csv(fnm)[,-1]
	params = strsplit(fnm, "-")[[1]][2:6]

	mse.mean = mean(output)
	mse.sd = sd(output)

	c(params,
	  mse=sprintf("%.2f", round(mse.mean, 2)),
	  sd=sprintf("%.2f", round(mse.sd, 2)))
})))


rownames(raw) = 1:nrow(raw)


names(raw) = c(param.names,
               "mse.mean",
               "mse.sd")

options(stringsAsFactors = FALSE)
raw = data.frame(apply(raw, 1:2, as.character))

raw = raw[!duplicated(raw[,1:5])] # can remove this line later
raw = dcast(setDT(raw), setup + n + p + sigma ~ alg, value.var=c("mse.mean", "mse.sd"))

raw = raw[order(as.numeric(raw$sigma)),]
raw = raw[order(as.numeric(raw$p)),]
raw = raw[order(as.numeric(raw$n)),]
raw = raw[order(as.numeric(raw$setup)),]
rownames(raw) = 1:nrow(raw)

write.csv(raw, file="output.csv")

#tab.all = cbind("", raw[,1 + c(2, 3, 5, 6, 10, 14, 7, 11, 15, 8, 12, 16, 9, 13, 17)])
#rmse.idx = c(5, 8, 11)
#for(iter in 1:nrow(tab.all)) {
#	best.idx = rmse.idx[which(as.numeric(tab.all[iter,rmse.idx]) == min(as.numeric(tab.all[iter,rmse.idx])))]
#	tab.all[iter,best.idx] = paste("\\bf", tab.all[iter,best.idx])
#}
#
#xtab.all = xtable(tab.all, align=c("r", "|", "|", "c", "|", rep("c", 3), "|", "|", rep("c", 3), "|", rep("c", 3), "|", rep("c", 3), "|", rep("c", 3)))
#print(xtab.all, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity,
#      hline.after = c(-1, -1, 0, 8, 16, 24, 32, 32), file = "simulation_results.tex")
