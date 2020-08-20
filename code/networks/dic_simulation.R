library(rjags)
library(parallel)
setwd("/Users/timvigers/Documents/GitHub/MS-Thesis")
set.seed(1017)

# Sample size
n = 160

# Structure 1: metabolite depends on T1D, methylation depends on T1D and metabolite
## T1D -> metab
alpha0 <- 5 # Intercept
alpha <- 2 # Slope
t1d <- rbinom(n,1,prob = 0.5) # Generate T1D cases and controls - our dataset is evenly split
metab <- rnorm(n, mean = alpha0 + alpha*t1d, sd = 1) # Outcome
## methyl | T1D, metab
int <- 10 # True intercept
beta <- 8 # T1D effect
gamma <- 5 # Metab effect

methyl <- rnorm(n, mean = int + beta*t1d + gamma*metab, sd = 1) # Outcome

df <- data.frame(t1d = t1d,metab = metab,methyl = methyl)

summary(lm(methyl~t1d+metab,df)) # Sanity check

# Permutation and MCMC parameters
nsim = 100
n_adapt = 1000
iter = 10000
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma")

# Permutation tests
dfs <- lapply(1:nsim,function(x){
  if (x == 1){
    d <- df
  } else {
    set.seed(x)
    d <- cbind.data.frame(df$t1d,sample(df$methyl),sample(df$metab))
  }
  colnames(d) <- colnames(df)
  return(d)
})
# Parallel - make cluster
no_cores = detectCores()
cl <- makeCluster(no_cores/4,type = "FORK")
# MCMC on each permuted set
mcmc_sim_perms <- parLapply(cl,dfs,function(x){
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
save(mcmc_sim_perms,file = "./data/networks/mcmc_sim_perms.Rdata")

# Try with scaled metabolites (like our real data)
df$metab <- scale(df$metab)
# Permutations
dfs <- lapply(1:nsim,function(x){
  if (x == 1){
    d <- df
  } else {
    set.seed(x)
    d <- cbind.data.frame(df$t1d,sample(df$methyl),sample(df$metab))
  }
  colnames(d) <- colnames(df)
  return(d)
})
# Parallel - make cluster
cl <- makeCluster(no_cores/4,type = "FORK")
# MCMC on each permuted set
mcmc_sim_perms_metab_scaled <- parLapply(cl,dfs,function(x){
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
save(mcmc_sim_perms_metab_scaled,file = "./data/networks/mcmc_sim_perms_metab_scaled.Rdata")

# Try with both scaled
df$methyl <- scale(df$methyl)
# Permutations
dfs <- lapply(1:nsim,function(x){
  if (x == 1){
    d <- df
  } else {
    set.seed(x)
    d <- cbind.data.frame(df$t1d,sample(df$methyl),sample(df$metab))
  }
  colnames(d) <- colnames(df)
  return(d)
})
# Parallel - make cluster
cl <- makeCluster(no_cores/4,type = "FORK")
# MCMC on each permuted set
mcmc_sim_perms_both_scaled <- parLapply(cl,dfs,function(x){
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
save(mcmc_sim_perms_both_scaled,file = "./data/networks/mcmc_sim_perms_both_scaled.Rdata")