library(rjags)
library(boot)
library(parallel)
# Load data
setwd("/Users/timvigers/GitHub/MS-Thesis")
load("./data/networks/pair_data.Rdata")
load("./data/networks/cits.Rdata")
# Bootstrap and MCMC parameters
n_adapt = 1000
iter = 10000
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma")
# Unique pairs from cit package
cits = cits[!(duplicated(cits[,c("methyl","metab")])),]
# Get DIC function
get_dic <- function(methyl,metab,data = pair_data,indices) {
  pair = data[,c("T1Dgroup",methyl,metab)]
  pair = pair[complete.cases(pair),]
  pair$T1Dgroup = ifelse(pair$T1Dgroup == "T1D control",0,1)
  d <- pair[indices,] 
  
  jags_data =
    list(t1d=c(d[,"T1Dgroup"]),methyl=c(d[,methyl]),
         metab=c(d[,metab]),N=nrow(d))
  
  dics =
    suppressWarnings(lapply(paste0("struct",1:24), function(x){
      mod =
        jags.model(paste0("./code/jags/",x,".jags"),quiet=T,
                   data = jags_data, n.adapt = n_adapt,n.chains = 2)
      dic =
        dic.samples(mod,n.iter = iter,progress.bar="none")
      return(round(sum(dic$deviance) + sum(dic$penalty),1))
    }))
  return(unlist(dics))
}
# Bootstrap
test <- boot(data = pair_data,statistic = get_dic,R=10,
             methyl = "cg01793374",metab = "gctof_235")
