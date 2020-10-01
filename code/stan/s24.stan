functions{
     int numLevels(int[] m) {
        int sorted[num_elements(m)];
        int count = 1;
        sorted = sort_asc(m);
        for (i in 2:num_elements(sorted)) {
          if (sorted[i] != sorted[i-1])
             count = count + 1;
        }
        return(count);
     }
}
data{
     // Define variables in data
     int<lower=1> Nobs;   // Number of observations (an integer)
     real methyl[Nobs];   // outcome variable
     real metab[Nobs];
     int<lower=0,upper=1> t1d[Nobs];   // outcome variable
}
transformed data{
     // Define transformed data
}
parameters{
     // Define parameters to estimate
     real<lower=0> sigma_methyl;
     real a_methyl;
     real b_metab_methyl;
     real<lower=0,upper=1> theta_t1d;
}
transformed parameters{
     // Transform parameters
     real mu_methyl[Nobs];
     for (i in 1:Nobs) {
        mu_methyl[i] = a_methyl + b_metab_methyl * metab[i];
     }
}
model{
     // Priors
     a_methyl ~ normal(0,100);
     b_metab_methyl ~ normal( 0, 10 );
     theta_t1d ~ beta(1, 1);

     // Likelihoods
     methyl ~ normal(mu_methyl, sigma_methyl);
     t1d ~ bernoulli(theta_t1d);
}
generated quantities {
     // simulate data from the posterior
     int<lower=0,upper=1> yrep_t1d[Nobs];
     // log-likelihood posterior
     vector[Nobs] log_lik_t1d;
     for (i in 1:num_elements(yrep_t1d)) {
       yrep_t1d[i] = binomial_rng(t1d[i], theta_t1d);
     }
     for (i in 1:Nobs) {
       log_lik_t1d[i] = binomial_lpmf(t1d[i] | 1, theta_t1d);
     }
}

