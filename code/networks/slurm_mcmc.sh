#!/bin/bash
#SBATCH -n 10
#SBATCH --mem=5000
#SBATCH --job-name=tim_mcmc
#SBATCH --output=tim_mcmc.log
Rscript $HOME/MS-Thesis/code/networks/mcmc.R