#!/bin/bash
#SBATCH -n 16
#SBATCH --job-name=tim_ewas_LI
#SBATCH --output=tim_ewas_LI.log
time Rscript $HOME/MS-Thesis/code/ewas_models_late_infancy.R
