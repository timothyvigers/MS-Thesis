library(rjags)
library(parallel)
setwd("/home/vigerst/MS-Thesis")

# Sample size
n = 160

# Permutation and MCMC parameters
nsim = 1000
n_adapt = 1000
iter = 10000
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma")

# True values
alpha0 <- 2 # Intercept
alpha <- 5 # Slope

# Multiple simulations
dfs <- lapply(1:nsim,function(x){
  set.seed(1016 + x)
  # Simulate data
  methyl <- rnorm(n)
  t1d <- rbinom(n,1,prob = plogis(alpha0 + alpha*methyl)) # Generate T1D cases and controls - our dataset is evenly split
  df <- data.frame(t1d = t1d,methyl = methyl)
  df$metab <- 0
  return(df)
})
# Parallel - make cluster
no_cores = detectCores()
cl <- makeCluster(no_cores/4,type = "FORK")
# MCMC on each permuted set
mcmc_simple_t1d <- parLapply(cl,dfs,function(x){
  jags_data =
    list(t1d=c(x$t1d),methyl=c(x$methyl),
         metab=c(x$metab),N=n)
  dics =
    lapply(paste0("struct",c(19,21,23)), function(x){
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
save(mcmc_simple_t1d,file = "./data/networks/mcmc_simple_t1d.Rdata")

# Two continuous
# Multiple simulations
dfs <- lapply(1:nsim,function(x){
  set.seed(1016 + x)
  # Simulate data
  methyl <- rnorm(n,5,2)
  metab <- rnorm(n,(alpha0 + alpha*methyl))
  df <- data.frame(metab = metab,methyl = methyl)
  df$t1d <- 0
  return(df)
})
# Parallel - make cluster
cl <- makeCluster(no_cores/4,type = "FORK")
# MCMC on each permuted set
mcmc_simple_omics <- parLapply(cl,dfs,function(x){
  jags_data =
    list(t1d=c(x$t1d),methyl=c(x$methyl),
         metab=c(x$metab),N=n)
  dics =
    lapply(paste0("struct",c(21,22,24)), function(x){
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
save(mcmc_simple_omics,file = "./data/networks/mcmc_simple_omics.Rdata")
# Scaled
dfs <- lapply(1:nsim,function(x){
  set.seed(1016 + x)
  # Simulate data
  methyl <- rnorm(n,5,2)
  metab <- rnorm(n,(alpha0 + alpha*methyl))
  df <- data.frame(metab = metab,methyl = methyl)
  df$t1d <- 0
  df$methyl <- scale(df$methyl)
  df$metab <- scale(df$metab)
  return(df)
})
# Parallel - make cluster
cl <- makeCluster(no_cores/4,type = "FORK")
# MCMC on each permuted set
mcmc_simple_omics_scaled <- parLapply(cl,dfs,function(x){
  jags_data =
    list(t1d=c(x$t1d),methyl=c(x$methyl),
         metab=c(x$metab),N=n)
  dics =
    lapply(paste0("struct",c(21,22,24)), function(x){
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
save(mcmc_simple_omics_scaled,file = "./data/networks/mcmc_simple_omics_scaled.Rdata")
