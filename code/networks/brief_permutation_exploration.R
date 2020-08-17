###########################
library(R2jags)

## set-up simple linear regression model and inits for jags
lm1_jags <- function(){
	# Likelihood:
	for (i in 1:N){
		b_mass[i] ~ dnorm(mu[i], tau) # tau is precision (1 / variance)
		mu[i] <- alpha + beta * b_length[i]
	}
	# Priors:
	alpha ~ dnorm(0, 0.01) # intercept
	beta ~ dnorm(0, 0.01) # slope
	sigma ~ dunif(0, 100) # standard deviation
	tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

init_values <- function(){
	list(alpha = rnorm(1), beta = rnorm(1), sigma = runif(1))
}

params <- c("alpha", "beta", "sigma")


###########################
### Slope not equal to 0 (=10)
set.seed(515) # Set a random seed for reproducibility of the simulation

samplesize <- 30 # Number of data points
b_length <- sort(rnorm(samplesize)) # Body length (explanatory variable)

int_true <- 30 # True intercept
slope_true <- 10 # True slope
mu <- int_true + slope_true * b_length # True means of normal distributions
sigma <- 5 # True standard deviation of normal distributions

b_mass <- rnorm(samplesize, mean = mu, sd = sigma) # Body mass (response variable)

snakes1 <- data.frame(b_length = b_length, b_mass = b_mass)
jagsdata_s1 <- with(snakes1, list(b_mass = b_mass, b_length = b_length, N = length(b_mass)))

fit_lm1 <- jags(data = jagsdata_s1, inits = init_values, parameters.to.save = params, model.file = lm1_jags,
			   n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = T)

summary(lm( b_mass ~ b_length, data=snakes1)) #check frequentist results

# Permute once
perm.vec <- NULL

for( i in 1:100 ){
	set.seed(515+i)
	snakes1p <- snakes1 #repeat data to permute
	snakes1p[,2] <- snakes1p[sample(1:30, replace=F),2]

	jagsdata_s1p <- with(snakes1p, list(b_mass = b_mass, b_length = b_length, N = length(b_mass)))

	lmres <- summary(lm( b_mass ~ b_length, data=snakes1p))

	fit_lm1p <- jags(data = jagsdata_s1p, inits = init_values, parameters.to.save = params, model.file = lm1_jags,
				   n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = T)
	perm.vec <- rbind(perm.vec, c(fit_lm1p$BUGSoutput$DIC, lmres$coef[2,4]) )
}

perm.vec #first column=DIC, second column=frequentist p-value
fit_lm1$BUGSoutput$DIC #compare perm.vec results to DIC for original data


###########################
### Slope equal to 0 (null case)
set.seed(6515) # Set a random seed for reproducibility of the simulation

samplesize <- 30 # Number of data points
b_length <- sort(rnorm(samplesize)) # Body length (explanatory variable)

int_true <- 30 # True intercept
slope_true <- 0 # True slope
mu <- int_true + slope_true * b_length # True means of normal distributions
sigma <- 5 # True standard deviation of normal distributions

b_mass <- rnorm(samplesize, mean = mu, sd = sigma) # Body mass (response variable)

snakes1 <- data.frame(b_length = b_length, b_mass = b_mass)
jagsdata_s1 <- with(snakes1, list(b_mass = b_mass, b_length = b_length, N = length(b_mass)))

fit_lm1_null <- jags(data = jagsdata_s1, inits = init_values, parameters.to.save = params, model.file = lm1_jags,
			   n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = T)

summary(lm( b_mass ~ b_length, data=snakes1)) #check frequentist results

# Permute once
perm.vec_null <- NULL

for( i in 1:100 ){
	set.seed(6515+i)
	snakes1p <- snakes1 #repeat data to permute
	snakes1p[,2] <- snakes1p[sample(1:30, replace=F),2]

	jagsdata_s1p <- with(snakes1p, list(b_mass = b_mass, b_length = b_length, N = length(b_mass)))

	lmres <- summary(lm( b_mass ~ b_length, data=snakes1p))

	fit_lm1p <- jags(data = jagsdata_s1p, inits = init_values, parameters.to.save = params, model.file = lm1_jags,
				   n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = T)
	perm.vec_null <- rbind(perm.vec_null, c(fit_lm1p$BUGSoutput$DIC, lmres$coef[2,4]) )
}

perm.vec_null
fit_lm1_null$BUGSoutput$DIC


###########################
### Slope not equal to 0 (=20)
set.seed(515) # Set a random seed for reproducibility of the simulation

samplesize <- 30 # Number of data points
b_length <- sort(rnorm(samplesize)) # Body length (explanatory variable)

int_true <- 30 # True intercept
slope_true <- 20 # True slope
mu <- int_true + slope_true * b_length # True means of normal distributions
sigma <- 5 # True standard deviation of normal distributions

b_mass <- rnorm(samplesize, mean = mu, sd = sigma) # Body mass (response variable)

snakes1 <- data.frame(b_length = b_length, b_mass = b_mass)
jagsdata_s1 <- with(snakes1, list(b_mass = b_mass, b_length = b_length, N = length(b_mass)))

fit_lm1_s20 <- jags(data = jagsdata_s1, inits = init_values, parameters.to.save = params, model.file = lm1_jags,
			   n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = T)

summary(lm( b_mass ~ b_length, data=snakes1)) #check frequentist results

# Permute once
perm.vec_s20 <- NULL

for( i in 1:100 ){
	set.seed(515+i)
	snakes1p <- snakes1 #repeat data to permute
	snakes1p[,2] <- snakes1p[sample(1:30, replace=F),2]

	jagsdata_s1p <- with(snakes1p, list(b_mass = b_mass, b_length = b_length, N = length(b_mass)))

	lmres <- summary(lm( b_mass ~ b_length, data=snakes1p))

	fit_lm1p <- jags(data = jagsdata_s1p, inits = init_values, parameters.to.save = params, model.file = lm1_jags,
				   n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = T)
	perm.vec_s20 <- rbind(perm.vec_s20, c(fit_lm1p$BUGSoutput$DIC, lmres$coef[2,4]) )
}

perm.vec_s20
fit_lm1_s20$BUGSoutput$DIC


###########################
### Slope = 10, resimulating each time
set.seed(515) # Set a random seed for reproducibility of the simulation

samplesize <- 30 # Number of data points
b_length <- sort(rnorm(samplesize)) # Body length (explanatory variable)

int_true <- 30 # True intercept
slope_true <- 10 # True slope
mu <- int_true + slope_true * b_length # True means of normal distributions
sigma <- 5 # True standard deviation of normal distributions

b_mass <- rnorm(samplesize, mean = mu, sd = sigma) # Body mass (response variable)

snakes1 <- data.frame(b_length = b_length, b_mass = b_mass)
jagsdata_s1 <- with(snakes1, list(b_mass = b_mass, b_length = b_length, N = length(b_mass)))

fit_lm1_resim <- jags(data = jagsdata_s1, inits = init_values, parameters.to.save = params, model.file = lm1_jags,
			   n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = T)

summary(lm( b_mass ~ b_length, data=snakes1))  

# Permute once
perm.vec_resim <- NULL

for( i in 1:100 ){
	set.seed(515+i)

	samplesize <- 30 # Number of data points
	b_length <- sort(rnorm(samplesize)) # Body length (explanatory variable)

	int_true <- 30 # True intercept
	slope_true <- 10 # True slope
	mu <- int_true + slope_true * b_length # True means of normal distributions
	sigma <- 5 # True standard deviation of normal distributions

	b_mass <- rnorm(samplesize, mean = mu, sd = sigma) # Body mass (response variable)
	snakes1p <- data.frame(b_length = b_length, b_mass = b_mass)

	jagsdata_s1p <- with(snakes1p, list(b_mass = b_mass, b_length = b_length, N = length(b_mass)))

	lmres <- summary(lm( b_mass ~ b_length, data=snakes1p))

	fit_lm1p <- jags(data = jagsdata_s1p, inits = init_values, parameters.to.save = params, model.file = lm1_jags,
				   n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = T)
	perm.vec_resim <- rbind(perm.vec_resim, c(fit_lm1p$BUGSoutput$DIC, lmres$coef[2,4]) )
}

perm.vec_resim
fit_lm1_resim$BUGSoutput$DIC

