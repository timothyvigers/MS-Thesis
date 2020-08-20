#!/bin/bash
#SBATCH -n 20
#SBATCH --job-name=tim_mcmc_unrelated_scaled
#SBATCH --output=tim_mcmc_unrelated_scaled.log
Rscript $HOME/MS-Thesis/code/networks/mcmc_permutations_random_pairs.R
