#!/bin/bash
#SBATCH --job-name=tim_datasets
#SBATCH --output=longitudinal_mediation_dataset.log
Rscript $HOME/MS-Thesis/code/mediation/make_longitudinal_mediation_dataset.R