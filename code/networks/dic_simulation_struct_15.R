library(rjags)
library(parallel)
setwd("/home/vigerst/MS-Thesis")

# Sample size
n = 160

# Permutation and MCMC parameters
nsim = 100
n_adapt = 1000
iter = 10000
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma")

# Structure 10: metabolite affect t1d and methylation
## True values
## T1D | metab
alpha0 <- 1 # Intercept
alpha <- 5 # Slope

## methyl | metab
gamma0 <- 10 # True intercept
gamma <- 5 # Metab effect

# Try with scaled metabolites (like our real data)
dfs <- lapply(1:nsim,function(x){
  set.seed(1016 + x)
  # Simulate data
  metab <- rnorm(n)
  t1d <- rbinom(n,1,prob = plogis(alpha0 + alpha*metab)) # Generate T1D cases and controls
  methyl <- rnorm(n, mean = gamma0 + gamma*metab, sd = 1)
  df <- data.frame(t1d = t1d,metab = metab,methyl = methyl)
  return(df)
})
# Parallel - make cluster
no_cores <- detectCores()
cl <- makeCluster(no_cores/4,type = "FORK")
# MCMC on each permuted set
mcmc_sim_struct15 <- parLapply(cl,dfs,function(x){
  jags_data =
    list(t1d=c(x$t1d),methyl=c(x$methyl),
         metab=c(x$metab),N=n)
  dics =
    lapply(paste0("struct",1:24), function(x){
      mod =
        jags.model(paste0("./code/jags/",x,".jags"),quiet=T,
                   data = jags_data, n.adapt = n_adapt,n.chains = 2)
      dic =
        dic.samples(mod,n.iter = iter,progress.bar="none")
      keep <- as.numeric(which(!is.nan(dic$penalty)))
      return(round(sum(dic$deviance[keep]) + sum(dic$penalty[keep]),1))
    })
  unlist(dics)
})
stopCluster(cl)
save(mcmc_sim_struct15,file = "./data/networks/mcmc_sim_struct15.Rdata")

# Try with both scaled
dfs <- lapply(1:nsim,function(x){
  set.seed(1016 + x)
  # Simulate data
  metab <- rnorm(n)
  t1d <- rbinom(n,1,prob = plogis(alpha0 + alpha*metab)) # Generate T1D cases and controls
  methyl <- rnorm(n, mean = gamma0 + gamma*metab, sd = 1)
  df <- data.frame(t1d = t1d,metab = metab,methyl = methyl)
  df$methyl <- scale(df$methyl)
  return(df)
})
# Parallel - make cluster
cl <- makeCluster(no_cores/4,type = "FORK")
# MCMC on each permuted set
mcmc_sim_struct15_both_scaled <- parLapply(cl,dfs,function(x){
  jags_data =
    list(t1d=c(x$t1d),methyl=c(x$methyl),
         metab=c(x$metab),N=n)
  dics =
    lapply(paste0("struct",1:24), function(x){
      mod =
        jags.model(paste0("./code/jags/",x,".jags"),quiet=T,
                   data = jags_data, n.adapt = n_adapt,n.chains = 2)
      dic =
        dic.samples(mod,n.iter = iter,progress.bar="none")
      keep <- as.numeric(which(!is.nan(dic$penalty)))
      return(round(sum(dic$deviance) + sum(dic$penalty),1))
    })
  unlist(dics)
})
stopCluster(cl)
save(mcmc_sim_struct15_both_scaled,file = "./data/networks/mcmc_sim_struct15_both_scaled.Rdata")