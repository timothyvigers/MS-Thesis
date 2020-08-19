#!/bin/bash
#SBATCH --job-name=tim_scaled_datasets
#SBATCH --output=tim_scaled_datasets.log
Rscript $HOME/MS-Thesis/code/networks/make_all_dataset.R
