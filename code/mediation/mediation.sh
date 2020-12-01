#!/bin/bash
#SBATCH -n 20
#SBATCH --job-name=tim_cits_mediation
#SBATCH --output=tim_cits_mediation.log
Rscript $HOME/MS-Thesis/code/mediation/mediation.R
