# Jags structures for R2jags
#1
struct1_jags <- function(){
  for (i in 1:N) {
    # metab | T1D
    metab[i] ~ dnorm(alpha0 + alpha*t1d[i], tau)
    # methyl | T1D, metab
    methyl[i] ~ dnorm(int + beta*t1d[i] + gamma*metab[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(methyl[i],int + beta*t1d[i] + gamma*metab[i], tau))
  }
  # Very flat priors
  alpha0 ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  int ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 2
struct2_jags <- function(){
  for (i in 1:N) {
    # methyl | T1D
    methyl[i] ~ dnorm(beta0 + beta*t1d[i], tau)
    # metab | T1D, methyl
    metab[i] ~ dnorm(int + alpha*t1d[i] + gamma*methyl[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(metab[i],int + alpha*t1d[i] + gamma*methyl[i], tau))
  }
  # Very flat priors
  alpha ~ dnorm(0,0.001)
  beta0 ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  int ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 3
struct3_jags <- function(){
  for (i in 1:N) {
    # metab | methyl
    metab[i] ~ dnorm(gamma0 + gamma*methyl[i], tau)
    # T1D | methyl, metab
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- int + alpha*metab[i] + beta*methyl[i]
    # Log likelihood for the network
    #LogLik[i] <- log(dbern(t1d[i],logit(int + alpha*t1d[i] + gamma*methyl[i]), tau))
  }
  # Very flat priors
  int ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  gamma0 ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}