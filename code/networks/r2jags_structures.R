# Jags structures for R2jags
struct1_jags <- function(){
  for (i in 1:N) {
    # metab | T1D
    metab[i] ~ dnorm(alpha0 + alpha*t1d[i], tau)
    # methyl | T1D, metab
    methyl[i] ~ dnorm(int + beta*t1d[i] + gamma*metab[i], tau)
    # Log likelihood
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
