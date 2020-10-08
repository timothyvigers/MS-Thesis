library(R2jags)
library(loo)
# Load data
setwd("/home/vigerst/MS-Thesis")
load("./data/networks/pair_data.Rdata")
load("./data/networks/cits.Rdata")
source("./code/networks/r2jags_structures.R")
# Permutation and MCMC parameters
cores = 24
n = 160
nsim = 100
vars = c("alpha0","alpha","beta0","beta","gamma0","gamma","LogLik")
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
# WAIC for each model
all_waics = lapply(dfs,function(x){
  jags_data =
    list(t1d=c(x$t1d),methyl=c(x$methyl),
         metab=c(x$metab),N=n)
  # R2jags - try all structures
  waics <- lapply(paste0("struct",1:24), function(w){
    j <- get(paste0(w,"_jags"))
    fit_lm <- 
      jags.parallel(data = jags_data,n.cluster = cores,
                            parameters.to.save = vars,model.file = j,
                            n.chains = 2, n.iter = 10000, n.burnin = 1000)
    loglik <- fit_lm$BUGSoutput$sims.list$LogLik
    waic <- suppressWarnings(waic(loglik))
    loo <- suppressWarnings(loo(loglik,cores = cores))
    return(c(waic$estimates["waic",1],loo$estimates["looic",1]))
  })
  waics <- as.data.frame(do.call(rbind,waics))
  colnames(waics) <- c("waic","looic")
  return(waics)
})
save(all_waics,file = "./data/networks/all_waic.Rdata")
