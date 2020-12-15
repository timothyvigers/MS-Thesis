#!/bin/bash
#SBATCH -n 24
#SBATCH --job-name=tim_cits_mediation
#SBATCH --output=tim_cits_mediation.log
time Rscript $HOME/MS-Thesis/code/mediation/cits_mediation_quasi.R
