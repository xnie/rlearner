rm(list = ls())

# the full dataset is available from
# https://github.com/gsbDBI/ExperimentData/tree/master/Mobilization/ProcessedData
data = read.csv("mobilization_no_unlisted 2.csv")

# W is intent to treat
# contact is received treatment

covariates = c("persons", "state", "county", "competiv",
               "st_sen", "st_hse", "newreg", "vote98",
               "vote00", "age", "female")
X = data[,which(names(data) %in% covariates)]
W = data$W
received_treatment = data$contact
Y = data$vote02

DF = data.frame(X, Y, W)

DF.nona = DF[!is.na(rowSums(DF)),]

idx.all = sample(c(sample(which(DF.nona$W == 0), sum(DF.nona$W) * 3/2), which(DF.nona$W == 1)))
DF.subset = DF.nona[idx.all,]

write.csv(DF.subset, "data_clean.csv", row.names = FALSE)
