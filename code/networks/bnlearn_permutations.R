library(bnlearn)
set.seed(1017)
# Load data
setwd("/home/vigerst/MS-Thesis")
load("./data/networks/pair_data.Rdata")
load("./data/networks/pair_list.Rdata")
load("./data/networks/cits.Rdata")
cits = cits[!(duplicated(cits[,c("methyl","metab")])),]
# Learn each model with T1D "continuous"
nperm <- 1000
bnlearn_perms = apply(cits,1,function(x){
  # Make dataset
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  pair = pair_data[,c("T1Dgroup",methyl,metab)]
  pair = pair[complete.cases(pair),]
  pair$T1Dgroup[pair$T1Dgroup == "T1D control"] = 0
  pair$T1Dgroup[pair$T1Dgroup == "T1D case"] = 1
  pair$T1Dgroup = as.numeric(pair$T1Dgroup) + rnorm(nrow(pair),sd = 0.0001)
  true_score <- score(hc(pair),pair)
  true_struct <- modelstring(hc(pair))
  perm_scores <- lapply(1:nperm, function(x){
    set.seed(1016+x)
    pair_perm <- pair
    pair_perm$T1Dgroup <- sample(pair_perm$T1Dgroup,replace = T)
    pair_perm[,methyl] <- sample(pair_perm[,methyl],replace = T)
    perm_score <- score(model2network(true_struct),pair_perm)
    perm_score
  })
  return(list(as.data.frame(unlist(perm_scores)),true_score))
})
save(bnlearn_perms,file = "./data/networks/bnlearn_perms.Rdata")
