library(rjags)
library(parallel)
set.seed(1017)
setwd("/home/vigerst/MS-Thesis")
# Load data
load("./data/networks/all_data.Rdata")
# Methyl and metabolite associations
gctof_step1 = read.csv("./data/candidate_selection/step_1/sv/gctof_SV_unadj_scaled.csv",stringsAsFactors = F)
hilic_step1 = read.csv("./data/candidate_selection/step_1/sv/hilic_SV_unadj_scaled.csv",stringsAsFactors = F)
lipid_step1 = read.csv("./data/candidate_selection/step_1/sv/lipid_SV_unadj_scaled.csv",stringsAsFactors = F)
oxylipin_step1 = read.csv("./data/candidate_selection/step_1/sv/oxylipin_SV_unadj_scaled.csv",stringsAsFactors = F)
vitd_step1 = read.csv("./data/candidate_selection/step_1/sv/vitd_SV_unadj_scaled.csv",stringsAsFactors = F)
step1 = rbind(gctof_step1,hilic_step1)
step1 = rbind(step1,lipid_step1)
step1 = rbind(step1,oxylipin_step1)
step1 = rbind(step1,vitd_step1)
step1 = metab_step3[order(metab_step3$p.value),]
# Import methylation association with T1D
methyl_step2 = read.csv("./data/candidate_selection/step_2/methyl_adj.csv",stringsAsFactors = F)
# Import metabolite associations with T1D
gctof_step3 = read.csv("./data/candidate_selection/step_3/gctof_adj.csv",stringsAsFactors = F)
hilic_step3 = read.csv("./data/candidate_selection/step_3/hilic_adj.csv",stringsAsFactors = F)
lipid_step3 = read.csv("./data/candidate_selection/step_3/lipid_adj.csv",stringsAsFactors = F)
oxylipin_step3 = read.csv("./data/candidate_selection/step_3/oxylipin_adj.csv",stringsAsFactors = F)
vitd_step3 = read.csv("./data/candidate_selection/step_3/vitd_adj.csv",stringsAsFactors = F)
metab_step3 = rbind(gctof_step3,hilic_step3)
metab_step3 = rbind(metab_step3,lipid_step3)
metab_step3 = rbind(metab_step3,oxylipin_step3)
metab_step3 = rbind(metab_step3,vitd_step3)
metab_step3 = metab_step3[order(metab_step3$p.value),]
# Get totally unrelated pairs
pairs = step1[step1$p.value >= 0.5,]
pairs = pairs[pairs$methyl %in% methyl_step2$methyl[methyl_step2$p.value >= 0.5],]
pairs = pairs[pairs$metab %in% metab_step3$metab[metab_step3$p.value >= 0.5],]
subset_n = 139
pairs = pairs[sample(1:nrow(pairs),subset_n),]
# DIC for each model with permutation tests 
# Permutation and MCMC parameters
nsim = 1
n_adapt = 1000
iter = 10000
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma")
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