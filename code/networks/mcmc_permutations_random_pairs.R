library(rjags)
library(parallel)
# Load data
setwd("/home/vigerst/MS-Thesis")
load("./data/networks/pair_data.Rdata")
load("./data/networks/pair_list.Rdata")
load("./data/networks/cits.Rdata")
load("./data/networks/all_data.Rdata")
set.seed(1017)
# Permutation and MCMC parameters
nsim = 100
n_adapt = 1000
iter = 10000
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma")
# DIC for each model with permutation tests 
# Get methyl and metabs that weren't in the pairs we selected
methyls <- colnames(all_data)[(!(colnames(all_data) %in% pairs$methyl))]
methyls <- methyls[grep("cg.*",methyls)]
metabs <- colnames(all_data)[!(colnames(all_data) %in% pairs$metab)]
metabs <- metabs[18:length(metabs)]
metabs <- metabs[grep("cg.*",metabs,invert = T)]
metabs <- metabs[grep("ch.*",metabs,invert = T)]
# Get data
subset_n = 139
pairs = pairs[1:subset_n,]
pairs$methyl <- sample(methyls,subset_n)
pairs$metab <- sample(metabs,subset_n)
# Iterate
all_perms = apply(pairs,1,function(x){
  methyl = as.character(x["methyl"])
  metab = as.character(x["metab"])
  print(paste(methyl,"and",metab))
  pair = all_data[,c("T1Dgroup",methyl,metab)]
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
  cl <- makeCluster(no_cores/4,type = "FORK")
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
subset_mcmc_perms = all_perms
save(subset_mcmc_perms,file = "./data/networks/subset_mcmc_perms.Rdata")