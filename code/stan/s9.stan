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
     int<lower=0,upper=1> t1d[Nobs];
}
transformed data{
     // Define transformed data
}
parameters{
     // Define parameters to estimate
     real<lower=0> sigma_methyl;
     real a_methyl;
     real b_t1d_methyl;
     real a_t1d;
     real b_metab_t1d;
}
transformed parameters{
     // Transform parameters
     real mu_methyl[Nobs];
     vector[Nobs] t1d;
     for(k in 1:Nobs) {
        t1d[k] = a_t1d + b_metab_t1d*metab[k];
     }
     for (i in 1:Nobs) {
        mu_methyl[i] = a_methyl + b_t1d_methyl * t1d[i];
     }
}
model{
     // Priors
     a_methyl ~ normal(0,100);
     b_t1d_methyl ~ normal( 0, 10 );
     a_t1d ~ normal( 0, 10 );
     b_metab_t1d ~ normal( 0, 10 );

     // Likelihoods
     methyl ~ normal(mu_methyl, sigma_methyl);
}
generated quantities {
     // simulate data from the posterior
     real yrep_methyl[Nobs];
     // log-likelihood posterior
     vector[Nobs] log_lik_methyl;
     for (i in 1:num_elements(yrep_methyl)) {
       yrep_methyl[i] = normal_rng(mu_methyl[i], sigma_methyl);
     }
     for (i in 1:Nobs) {
       log_lik_methyl[i] = normal_lpdf(methyl[i] | mu_methyl[i], sigma_methyl);
     }
}

