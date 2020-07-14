library(rjags)
library(parallel)
# Load data
setwd("/home/vigerst/MS-Thesis")
load("./data/networks/pair_data.Rdata")
load("./data/networks/cits.Rdata")
# Permutation and MCMC parameters
nsim = 100
n_adapt = 1000
iter = 10000
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma")
# Unique pairs from cit package
cits = cits[!(duplicated(cits[,c("methyl","metab")])),]
# DIC for each model with permutation tests
all_perms = apply(cits,1,function(x){
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  pair = pair_data[,c("T1Dgroup",methyl,metab)]
  pair = pair[complete.cases(pair),]
  pair$T1Dgroup = ifelse(pair$T1Dgroup == "T1D control",0,1)
  # Permutation datasets
  dfs <- lapply(1:nsim,function(x){
    if (x == 1){
      df <- pair
    } else {
      set.seed(x)
      df <- cbind.data.frame(pair$T1Dgroup,sample(pair[,methyl]),sample(pair[,metab]))
    }
    colnames(df) <- colnames(pair)
    return(df)
  })
  # Parallel - make cluster
  no_cores = detectCores()
  cl <- makeCluster(no_cores,type = "FORK")
  # MCMC on each permuted set
  mcmc_perms <- parLapply(cl=cl,dfs,function(x){
    jags_data =
      list(t1d=c(x[,"T1Dgroup"]),methyl=c(x[,methyl]),
           metab=c(x[,metab]),N=nrow(pair))
    dics =
      suppressWarnings(lapply(paste0("struct",1:24), function(x){
        mod =
          jags.model(paste0("./code/jags/",x,".jags"),quiet=T,
                     data = jags_data, n.adapt = n_adapt,n.chains = 2)
        dic =
          dic.samples(mod,n.iter = iter,progress.bar="none")
        return(round(sum(dic$deviance) + sum(dic$penalty),1))
      }))
    unlist(dics)
  })
  stopCluster(cl)
  return(mcmc_perms)
})
save(all_perms,file = "./data/networks/all_mcmc_perms.Rdata")