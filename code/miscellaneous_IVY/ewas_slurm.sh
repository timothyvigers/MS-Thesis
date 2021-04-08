#!/bin/bash
#SBATCH -n 16
#SBATCH --job-name=tim_ewas
#SBATCH --output=tim_ewas.log
time Rscript $HOME/MS-Thesis/code/ewas_models_late_infancy.R
time Rscript $HOME/MS-Thesis/code/ewas_models_childhood.R