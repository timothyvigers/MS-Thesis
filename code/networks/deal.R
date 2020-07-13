library(deal)
# Load data
setwd("/Users/timvigers/GitHub/MS-Thesis/")
load("./data/networks/pair_data.Rdata")
load("./data/networks/cits.Rdata")
# Unique pairs from cit package
cits = cits[!(duplicated(cits[,c("methyl","metab")])),]
# Iterate through pairs
# Testing
x <- cits[1,]
methyl = as.character(x["methyl"])
metab = as.character(x["metab"])
pair = pair_data[,c("T1Dgroup",methyl,metab)]
pair = pair[complete.cases(pair),]
pair$T1Dgroup <- factor(pair$T1Dgroup,levels = c("T1D control","T1D case"))
# Deal package
net <- network(pair,specifygraph=FALSE,inspectprob=FALSE,doprob=TRUE)
l <- learn(net,pair)
