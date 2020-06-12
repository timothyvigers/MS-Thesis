library(rjags)
# Load data
setwd("/Users/timvigers/GitHub/MS-Thesis")
load("./data/networks/pair_data.Rdata")
load("./data/networks/cits.Rdata")
# MCMC parameters
n_adapt = 1000
iter = 10000
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma")
# Permutations
nsim = 1000
# All pairs from cit package
# DIC for each model
all_perm_dics = apply(cits[1:2,],1,function(x){
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  temp = pair_data[,c("T1Dgroup",methyl,metab)]
  temp = temp[complete.cases(temp),]
  temp$T1Dgroup = ifelse(temp$T1Dgroup == "T1D control",0,1)
  N = nrow(temp)
  # Permutation test, return the minimum each time
  perm_structs = 
    lapply(1:nsim, function(x){
      perm = sample(nrow(temp))
      jags_data = 
        list(t1d=c(temp[perm,"T1Dgroup"]),methyl=c(temp[,methyl]),
             metab=c(temp[,metab]),
             N=N)
      dics = 
        suppressWarnings(lapply(paste0("struct",1:24), function(x){
          mod = 
            jags.model(paste0("./code/jags/",x,".jags"),quiet=T,
                       data = jags_data, n.adapt = n_adapt,n.chains = 2)
          dic = 
            dic.samples(mod,n.iter = iter,progress.bar="none")
          return(round(sum(dic$deviance) + sum(dic$penalty),1))
        }))
      which.min(dics)
    })
save(all_perm_dics,file = "./data/networks/all_perm_dic.Rdata")