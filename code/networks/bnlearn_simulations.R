library(bnlearn)
set.seed(1017)
setwd("C:/Users/timbv/Documents/GitHub/MS-Thesis")
## Simulate structure 1:
nsim = 10000
n = 160
correct <- "[t1d][metab|t1d][methyl|t1d:metab]"
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
all_learned = lapply(dfs,function(x){
  # Add noise to T1D
  x$t1d = as.numeric(x$t1d) + rnorm(n,sd = 0.0001)
  learned <- hc(x)
  mod <- modelstring(learned)
  check <- ifelse(mod == correct,"Correct",
                  abs(score(model2network(correct),x)-score(learned,x)))
  return(check)
})
struct1_learned_sim <- as.data.frame(unlist(all_learned))
save(struct1_learned_sim,file = "./data/networks/struct1_learned_sim.Rdata")
# Structure 15: metabolite affect t1d and methylation
correct <- "[metab][t1d|metab][methyl|metab]"
## True values
## T1D | metab
alpha0 <- 5 # Intercept
alpha <- 2 # Slope
## methyl | metab
gamma0 <- 10 # True intercept
gamma <- 8 # Metab effect
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
all_learned = lapply(dfs,function(x){
  # Add noise to T1D
  x$t1d = as.numeric(x$t1d) + rnorm(n,sd = 0.0001)
  learned <- hc(x)
  mod <- modelstring(learned)
  check <- ifelse(mod == correct,"Correct",
                  abs(score(model2network(correct),x)-score(learned,x)))
  return(check)
})
struct10_learned_sim <- as.data.frame(unlist(all_learned))
save(struct10_learned_sim,file = "./data/networks/struct10_learned_sim.Rdata")

