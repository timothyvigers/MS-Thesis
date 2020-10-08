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
     real methyl[Nobs];
     real metab[Nobs];   // outcome variable
     int<lower=0,upper=1> t1d[Nobs];   // outcome variable
}
transformed data{
     // Define transformed data
}
parameters{
     // Define parameters to estimate
     real<lower=0> sigma_metab;
     real a_metab;
     real b_methyl_metab;
     real a_t1d;
     real b_methyl_t1d;
}
transformed parameters{
     // Transform parameters
     real mu_metab[Nobs];
     real<lower=0,upper=1> theta_t1d[Nobs];
     for (i in 1:Nobs) {
        theta_t1d[i] = a_t1d + b_methyl_t1d * methyl[i];
     }
     for (i in 1:Nobs) {
        mu_metab[i] = a_metab + b_methyl_metab * methyl[i];
     }
}
model{
     // Priors
     a_metab ~ normal(0,100);
     b_methyl_metab ~ normal( 0, 10 );
     a_t1d ~ normal(0,100);
     b_methyl_t1d ~ normal( 0, 10 );

     // Likelihoods
     metab ~ normal(mu_metab, sigma_metab);
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

