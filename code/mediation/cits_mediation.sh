#!/bin/bash
#SBATCH -n 12
#SBATCH --job-name=tim_cits_mediation
#SBATCH --output=tim_cits_mediation.log
Rscript $HOME/MS-Thesis/code/mediation/cits_mediation.R
