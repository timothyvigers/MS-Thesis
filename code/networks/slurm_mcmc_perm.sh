#!/bin/bash
#SBATCH -n 20
#SBATCH --job-name=tim_mcmc
#SBATCH --output=tim_mcmc.log
Rscript $HOME/MS-Thesis/code/networks/mcmc_permutations_random_pairs.R
