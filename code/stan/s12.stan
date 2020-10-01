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
     int<lower=0,upper=1> t1d[Nobs];
}
transformed data{
     // Define transformed data
}
parameters{
     // Define parameters to estimate
     real a_methyl;
     real b_t1d_methyl;
     real<lower=0> sigma_metab;
     real a_metab;
     real b_methyl_metab;
}
transformed parameters{
     // Transform parameters
     vector[Nobs] methyl;
     real mu_metab[Nobs];
     for(k in 1:Nobs) {
        methyl[k] = a_methyl + b_t1d_methyl*t1d[k];
     }
     for (i in 1:Nobs) {
        mu_metab[i] = a_metab + b_methyl_metab * methyl[i];
     }
}
model{
     // Priors
     a_methyl ~ normal( 0, 10 );
     b_t1d_methyl ~ normal( 0, 10 );
     a_metab ~ normal(0,100);
     b_methyl_metab ~ normal( 0, 10 );

     // Likelihoods
     metab ~ normal(mu_metab, sigma_metab);
}
generated quantities {
     // simulate data from the posterior
     real yrep_metab[Nobs];
     // log-likelihood posterior
     vector[Nobs] log_lik_metab;
     for (i in 1:num_elements(yrep_metab)) {
       yrep_metab[i] = normal_rng(mu_metab[i], sigma_metab);
     }
     for (i in 1:Nobs) {
       log_lik_metab[i] = normal_lpdf(metab[i] | mu_metab[i], sigma_metab);
     }
}

