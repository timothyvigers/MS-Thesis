library(rjags)
library(parallel)
setwd("/home/vigerst/MS-Thesis")
set.seed(1017)

# Sample size
n = 160

# Permutation and MCMC parameters
nsim = 1000
n_adapt = 1000
iter = 10000
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma")

# Structure 1: metabolite depends on T1D, methylation depends on T1D and metabolite
## True values
## T1D -> metab
alpha0 <- 5 # Intercept
alpha <- 2 # Slope

## methyl | T1D, metab
int <- 10 # True intercept
beta <- 8 # T1D effect
gamma <- 5 # Metab effect

# Multiple simulations
dfs <- lapply(1:nsim,function(x){
  set.seed(1016 + x)
  # Simulate data
  t1d <- rbinom(n,1,prob = 0.5) # Generate T1D cases and controls - our dataset is evenly split
  metab <- rnorm(n, mean = alpha0 + alpha*t1d, sd = 1) # Outcome
  methyl <- rnorm(n, mean = int + beta*t1d + gamma*metab, sd = 1) # Outcome
  df <- data.frame(t1d = t1d,metab = metab,methyl = methyl)
  return(df)
})
# Parallel - make cluster
no_cores = detectCores()
cl <- makeCluster(no_cores/4,type = "FORK")
# MCMC on each permuted set
mcmc_sim <- parLapply(cl,dfs,function(x){
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
      return(round(sum(dic$deviance) + sum(dic$penalty),1))
    })
  unlist(dics)
})
stopCluster(cl)
save(mcmc_sim,file = "./data/networks/mcmc_sim.Rdata")

# Try with scaled metabolites (like our real data)
# Permutations
dfs <- lapply(1:nsim,function(x){
  set.seed(1016 + x)
  # Simulate data
  t1d <- rbinom(n,1,prob = 0.5) # Generate T1D cases and controls - our dataset is evenly split
  metab <- rnorm(n, mean = alpha0 + alpha*t1d, sd = 1) # Outcome
  methyl <- rnorm(n, mean = int + beta*t1d + gamma*metab, sd = 1) # Outcome
  df <- data.frame(t1d = t1d,metab = metab,methyl = methyl)
  df$metab <- scale(df$metab)
  return(df)
})
# Parallel - make cluster
cl <- makeCluster(no_cores/4,type = "FORK")
# MCMC on each permuted set
mcmc_sim_metab_scaled <- parLapply(cl,dfs,function(x){
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
      return(round(sum(dic$deviance) + sum(dic$penalty),1))
    })
  unlist(dics)
})
stopCluster(cl)
save(mcmc_sim_metab_scaled,file = "./data/networks/mcmc_sim_metab_scaled.Rdata")

# Try with both scaled
dfs <- lapply(1:nsim,function(x){
  set.seed(1016 + x)
  # Simulate data
  t1d <- rbinom(n,1,prob = 0.5) # Generate T1D cases and controls - our dataset is evenly split
  metab <- rnorm(n, mean = alpha0 + alpha*t1d, sd = 1) # Outcome
  methyl <- rnorm(n, mean = int + beta*t1d + gamma*metab, sd = 1) # Outcome
  df <- data.frame(t1d = t1d,metab = metab,methyl = methyl)
  df$metab <- scale(df$metab)
  df$methyl <- scale(df$methyl)
  return(df)
})
# Parallel - make cluster
cl <- makeCluster(no_cores/4,type = "FORK")
# MCMC on each permuted set
mcmc_sim_both_scaled <- parLapply(cl,dfs,function(x){
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
      return(round(sum(dic$deviance) + sum(dic$penalty),1))
    })
  unlist(dics)
})
stopCluster(cl)
save(mcmc_sim_both_scaled,file = "./data/networks/mcmc_sim_both_scaled.Rdata")