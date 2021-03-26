#!/bin/bash
#SBATCH -n 16
#SBATCH --job-name=tim_mediation
#SBATCH --output=tim_mediation.log
Rscript $HOME/MS-Thesis/code/ewas_models_late_infancy.R
Rscript $HOME/MS-Thesis/code/ewas_models_childhood.R