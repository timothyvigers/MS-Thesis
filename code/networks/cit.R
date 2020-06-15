library(cit)
# Load list and data
setwd("/Users/timvigers/GitHub/MS-Thesis")
load("./data/networks/pair_list.Rdata")
load("./data/networks/pair_data.Rdata")
pair_data$T1Dgroup = factor(pair_data$T1Dgroup)
# Run CIT package
cits = apply(pairs, 1, function(x){
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  temp = pair_data[,c("T1Dgroup",methyl,metab)]
  temp = temp[complete.cases(temp),]
  y = as.numeric(temp$T1Dgroup) - 1
  X = temp[,c(methyl,metab)]
  cinf1 = cit.bp(L = X[,1],G = X[,2],T = y)
  cinf2 = cit.bp(L = X[,2],G = X[,1],T = y)
  d = rbind(cinf1,cinf2)
  d$direction = c(">","<")
  d$methyl = methyl
  d$metab = metab
  return(d[,c("methyl","direction","metab","p_cit")])
})
# Convert to single dataframe, remove non-significant rows
cits = as.data.frame(do.call(rbind,cits))
cits = cits[cits$p_cit < 0.05,]
rownames(cits) = c(1:nrow(cits))
# Write
save(cits,file = "./data/networks/cits.Rdata")