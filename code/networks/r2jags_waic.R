library(R2jags)
library(parallel)
# Load data
setwd("C:/Users/timbv/Documents/GitHub/MS-Thesis")
load("./data/networks/pair_data.Rdata")
load("./data/networks/cits.Rdata")
set.seed(1017)
# Permutation and MCMC parameters
nsim = 100
n_adapt = 1000
iter = 10000
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma","LogLik")
# Unique pairs from cit package
cits = cits[!(duplicated(cits[,c("methyl","metab")])),]
# DIC for each model
all_waics = apply(cits[1:5,],1,function(x){
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  temp = pair_data[,c("T1Dgroup",methyl,metab)]
  temp = temp[complete.cases(temp),]
  temp$T1Dgroup = ifelse(temp$T1Dgroup == "T1D control",0,1)
  N = nrow(temp)
  jags_data = list(t1d=c(temp[,"T1Dgroup"]),methyl=c(temp[,methyl]),
                   metab=c(temp[,metab]),
                   N=N)
  # R2jags - try all structures
  waics <- lapply(paste0("struct",1:3), function(w){
    j <- get(paste0(w,"_jags"))
    fit_lm <- suppressWarnings(
      jags(data = jags_data,
           parameters.to.save = vars,model.file = j,
           n.chains = 2, n.iter = iter, n.burnin = n_adapt))
    loglik <- fit_lm$BUGSoutput$sims.list$LogLik
    waic <- suppressWarnings(waic(loglik))
    loo <- suppressWarnings(loo(loglik))
    return(c(waic$estimates["waic",1],loo$estimates["looic",1]))
  })
  return(waics)
})
save(all_waics,file = "./data/networks/all_waic.Rdata")
