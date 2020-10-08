library(bnlearn)
# Load data
setwd("C:/Users/timbv/Documents/GitHub/MS-Thesis")
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
  pair$T1Dgroup[pair$T1Dgroup == "T1D control"] = 0
  pair$T1Dgroup[pair$T1Dgroup == "T1D case"] = 1
  pair$T1Dgroup = as.numeric(pair$T1Dgroup) + rnorm(nrow(pair),sd = 0.0001)
  learned <- hc(pair)
  return(learned)
})
save(all_learned,file = "./data/networks/bnlearn_continuous_t1d.Rdata")
# Same again for cit-selected only
cit_learned = apply(cits,1,function(x){
  # Make dataset
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  pair = pair_data[,c("T1Dgroup",methyl,metab)]
  pair = pair[complete.cases(pair),]
  pair$T1Dgroup[pair$T1Dgroup == "T1D control"] = 0
  pair$T1Dgroup[pair$T1Dgroup == "T1D case"] = 1
  pair$T1Dgroup = as.numeric(pair$T1Dgroup) + rnorm(nrow(pair),sd = 0.0001)
  learned <- hc(pair)
  return(learned)
})
save(cit_learned,file = "./data/networks/bnlearn_cit_continuous_t1d.Rdata")
