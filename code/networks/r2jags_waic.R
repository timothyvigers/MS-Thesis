library(R2jags)
library(parallel)
# Load data
setwd("/Users/timvigers/Documents/GitHub/MS-Thesis")
load("./data/networks/pair_data.Rdata")
load("./data/networks/cits.Rdata")
set.seed(1017)
# MCMC parameters
n_adapt = 1000
iter = 10000
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma","LogLik")
# Unique pairs from cit package
cits = cits[!(duplicated(cits[,c("methyl","metab")])),]
# Parallel
# Make cluster
no_cores = 10
cl <- makeCluster(no_cores,type = "FORK")
# DIC for each model
all_dics = parApply(cl=cl,cits,1,function(x){
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  temp = pair_data[,c("T1Dgroup",methyl,metab)]
  temp = temp[complete.cases(temp),]
  temp$T1Dgroup = ifelse(temp$T1Dgroup == "T1D control",0,1)
  N = nrow(temp)
  jags_data = list(t1d=c(temp[,"T1Dgroup"]),methyl=c(temp[,methyl]),
                   metab=c(temp[,metab]),
                   N=N)
  # R2jags
  fit_lm1 <- jags(data = jags_data, parameters.to.save = vars, model.file = struct1_jags,
                  n.chains = 2, n.iter = 10000, n.burnin = 1000)
  # data for jags
 
  # Test all structures
  dics = lapply(paste0("struct",1:24), function(x){
    mod = jags.model(paste0("./code/jags/",x,".jags"),quiet=T,
                     data = jags_data, n.adapt = n_adapt,n.chains = 2)
    dic = dic.samples(mod,n.iter = iter,progress.bar = "none")
    return(round(sum(dic$deviance) + sum(dic$penalty),1))
  })
  return(unlist(dics))
})
save(all_dics,file = "./data/networks/all_dic.Rdata")
# Samples for all models
all_samples = parApply(cl=cl,cits,1,function(x){
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  temp = pair_data[,c("T1Dgroup",methyl,metab)]
  temp = temp[complete.cases(temp),]
  temp$T1Dgroup = ifelse(temp$T1Dgroup == "T1D control",0,1)
  # data for jags
  N = nrow(temp)
  jags_data = list(t1d=c(temp[,"T1Dgroup"]),methyl=c(temp[,methyl]),
                   metab=c(temp[,metab]),
                   N=N)
  # Test all structures
  samples = lapply(paste0("struct",1:24), function(x){
    mod = jags.model(paste0("./code/jags/",x,".jags"),quiet=T,
                     data = jags_data, n.adapt = n_adapt,n.chains = 2)
    coda = coda.samples(mod,variable.names = vars,n.iter = iter,progress.bar="none")
    return(coda)
  })
  return(samples)
})
save(all_samples,file = "./data/networks/all_samples.Rdata")