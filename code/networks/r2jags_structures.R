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
    LogLik[i] <- log(dbin(t1d[i],q[i],1))
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
# 4
struct4_jags <- function(){
  for (i in 1:N) {
    # T1D | methyl
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- beta0 + beta*methyl[i]
    # metab | T1D, methyl
    metab[i] ~ dnorm(int + alpha*t1d[i] + gamma*methyl[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(metab[i],int + alpha*t1d[i] + gamma*methyl[i],tau))
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
# 5
struct5_jags <- function(){
  for (i in 1:N) {
    # methyl | metab
    methyl[i] ~ dnorm(gamma0 + gamma*metab[i], tau)
    # T1D | methyl, metab
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- int + beta*methyl[i] + alpha*metab[i]
    # Log likelihood for the network
    LogLik[i] <- log(dbin(t1d[i],q[i],1))
  }
  # Very flat priors
  alpha ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  int ~ dnorm(0,0.001)
  gamma0 ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 6
struct6_jags <- function(){
  for (i in 1:N) {
    # T1D | metab
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- alpha0 + alpha*metab[i]
    # methyl | metab, T1D
    methyl[i] ~ dnorm(int + gamma*metab[i] + beta*t1d[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(methyl[i],int + gamma*metab[i] + beta*t1d[i], tau))
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
# 7
struct7_jags <- function(){
  for (i in 1:N) {
    # methyl | T1D
    methyl[i] ~ dnorm(beta0 + beta*t1d[i], tau)
    # metab | T1D
    metab[i] ~ dnorm(alpha0 + alpha*t1d[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(metab[i],alpha0 + alpha*t1d[i], tau)) + 
      log(dnorm(methyl[i],beta0 + beta*t1d[i], tau))
  }
  # Very flat priors
  alpha0 ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
  beta0 ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 8
struct8_jags <- function(){
  for (i in 1:N) {
    # T1D | methyl
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- beta0 + beta*methyl[i]
    # metab | T1D
    metab[i] ~ dnorm(alpha0 + alpha*t1d[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(metab[i],alpha0 + alpha*t1d[i], tau))
  }
  # Very flat priors
  alpha0 ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
  beta0 ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 9
struct9_jags <- function(){
  for (i in 1:N) {
    # T1D | metab
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- alpha0 + alpha*metab[i]
    # methyl | T1D
    methyl[i] ~ dnorm(beta0 + beta*t1d[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(methyl[i],beta0 + beta*t1d[i], tau))
  }
  # Very flat priors
  alpha0 ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
  beta0 ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 10
struct10_jags <- function(){
  for (i in 1:N) {
    # T1D | metab, methyl
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- int + alpha*metab[i] + beta*methyl[i]
    # Log likelihood for the network
    LogLik[i] <- log(dbin(t1d[i],q[i],1))
  }
  # Very flat priors
  int ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 11
struct11_jags <- function(){
  for (i in 1:N) {
    # T1D | methyl
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- beta0 + beta*methyl[i]
    # metab | methyl
    metab[i] ~ dnorm(gamma0 + gamma*methyl[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(metab[i],gamma0 + gamma*methyl[i], tau)) + 
      log(dbin(t1d[i],q[i],1))
  }
  # Very flat priors
  beta0 ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  gamma0 ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 12
struct12_jags <- function(){
  for (i in 1:N) {
    # methyl | T1D
    methyl[i] ~ dnorm(beta0 + beta*t1d[i], tau)
    # metab | methyl
    metab[i] ~ dnorm(gamma0 + gamma*methyl[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(metab[i],gamma0 + gamma*methyl[i], tau))
  }
  # Very flat priors
  beta0 ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  gamma0 ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 13
struct13_jags <- function(){
  for (i in 1:N) {
    # methyl | metab
    methyl[i] ~ dnorm(gamma0 + gamma*metab[i], tau)
    # T1D | methyl
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- beta0 + beta*methyl[i]
    # Log likelihood for the network
    LogLik[i] <- log(dbin(t1d[i],q[i],1))
  }
  # Very flat priors
  beta0 ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  gamma0 ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 14
struct14_jags <- function(){
  for (i in 1:N) {
    # methyl | metab, T1D
    methyl[i] ~ dnorm(int + gamma*metab[i] + beta*t1d[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(methyl[i],int + gamma*metab[i] + beta*t1d[i], tau))
  }
  # Very flat priors
  int ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 15
struct15_jags <- function(){
  for (i in 1:N) {
    # methyl | metab
    methyl[i] ~ dnorm(gamma0 + gamma*metab[i], tau)
    # T1D | metab
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- alpha0 + alpha*metab[i]
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(methyl[i],gamma0 + gamma*metab[i], tau)) +
      log(dbin(t1d[i],q[i],1))
  }
 # Very flat priors
  alpha0 ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
  gamma0 ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 16
struct16_jags <- function(){
  for (i in 1:N) {
    # metab | T1D
    metab[i] ~ dnorm(alpha0 + alpha*t1d[i], tau)
    # methyl | metab
    methyl[i] ~ dnorm(gamma0 + gamma*metab[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(methyl[i],gamma0 + gamma*metab[i], tau))
  }
  # Very flat priors
  alpha0 ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
  gamma0 ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 17
struct17_jags <- function(){
  for (i in 1:N) {
    # metab | methyl
    metab[i] ~ dnorm(gamma0 + gamma*methyl[i], tau)
    # T1D | metab
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- alpha0 + alpha*metab[i]
    # Log likelihood for the network
    LogLik[i] <- log(dbin(t1d[i],q[i],1))
  }
  # Very flat priors
  alpha0 ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
  gamma0 ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 18
struct18_jags <- function(){
  for (i in 1:N) {
    # metab | methyl, T1D
    metab[i] ~ dnorm(int + gamma*methyl[i] + alpha*t1d[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(metab[i],int + gamma*methyl[i] + alpha*t1d[i], tau))
  }
  # Very flat priors
  int ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 19
struct19_jags <- function(){
  for (i in 1:N) {
    # methyl | T1D
    methyl[i] ~ dnorm(beta0 + beta*t1d[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(methyl[i],beta0 + beta*t1d[i], tau))
  }
  # Very flat priors
  beta0 ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 20
struct20_jags <- function(){
  for (i in 1:N) {
    # metab | T1D
    metab[i] ~ dnorm(alpha0 + alpha*t1d[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(metab[i],alpha0 + alpha*t1d[i], tau))
  }
  # Very flat priors
  alpha0 ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 21
struct21_jags <- function(){
  for (i in 1:N) {
    # T1D | methyl
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- beta0 + beta*methyl[i]
    # Log likelihood for the network
    LogLik[i] <- log(dbin(t1d[i],q[i],1))
  }
  # Very flat priors
  beta0 ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
}
# 22
struct22_jags <- function(){
  for (i in 1:N) {
    # metab | methyl
    metab[i] ~ dnorm(gamma0 + gamma*methyl[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(metab[i],gamma0 + gamma*methyl[i], tau))
  }
  # Very flat priors
  gamma0 ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}
# 23
struct23_jags <- function(){
  for (i in 1:N) {
    # T1D | metab
    t1d[i] ~ dbern(q[i])
    logit(q[i]) <- alpha0 + alpha*metab[i]
    # Log likelihood for the network
    LogLik[i] <- log(dbin(t1d[i],q[i],1))
  }
  # Very flat priors
  alpha0 ~ dnorm(0,0.001)
  alpha ~ dnorm(0,0.001)
}
# 24
struct24_jags <- function(){
  # methyl | metab
  for (i in 1:N) {
    methyl[i] ~ dnorm(gamma0 + gamma*metab[i], tau)
    # Log likelihood for the network
    LogLik[i] <- log(dnorm(methyl[i],gamma0 + gamma*metab[i], tau))
  }
  # Very flat priors
  gamma0 ~ dnorm(0,0.001)
  gamma ~ dnorm(0,0.001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1/tau
}