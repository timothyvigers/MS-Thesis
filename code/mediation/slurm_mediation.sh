#!/bin/bash
#SBATCH -n 16
#SBATCH --job-name=tim_selection
#SBATCH --output=tim_selection.log
time Rscript $HOME/MS-Thesis/code/mediation/longitudinal_selection.R