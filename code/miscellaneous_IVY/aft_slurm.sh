#!/bin/bash
#SBATCH -n 16
#SBATCH --job-name=tim_aft
#SBATCH --output=tim_aft.log
time Rscript $HOME/MS-Thesis/code/miscellaneous_IVY/mediation_accelerated_failure.R