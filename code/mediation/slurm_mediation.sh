#!/bin/bash
#SBATCH -n 8
#SBATCH --job-name=tim_mediation
#SBATCH --output=tim_mediation.log
time Rscript $HOME/MS-Thesis/code/mediation/mediation_regmedint.R