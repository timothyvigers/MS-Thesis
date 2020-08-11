library(bnlearn)
# Load data
setwd("/Users/timvigers/Documents/GitHub/MS-Thesis/")
load("./data/networks/pair_data.Rdata")
load("./data/networks/pair_list.Rdata")
load("./data/networks/cits.Rdata")
cits = cits[!(duplicated(cits[,c("methyl","metab")])),]
set.seed(1017)
# Learn each model with continuous variables discretized
all_learned = apply(pairs,1,function(x){
  # Make dataset
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  pair = pair_data[,c("T1Dgroup",methyl,metab)]
  pair = pair[complete.cases(pair),]
  pair$T1Dgroup = factor(pair$T1Dgroup)
  d <- discretize(pair,breaks = 5,method = "interval")
  learned <- hc(d,score = "aic")
  return(learned)
})
save(all_learned,file = "./data/networks/bnlearn_discretized.Rdata")
# Same again for cit-selected only
cit_learned = apply(cits,1,function(x){
  # Make dataset
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  pair = pair_data[,c("T1Dgroup",methyl,metab)]
  pair = pair[complete.cases(pair),]
  pair$T1Dgroup = factor(pair$T1Dgroup)
  d <- discretize(pair,breaks = 5,method = "interval")
  learned <- hc(d,score = "aic")
  return(learned)
})
save(cit_learned,file = "./data/networks/bnlearn_cit_discretized.Rdata")