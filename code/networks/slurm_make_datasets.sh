#!/bin/bash
#SBATCH --job-name=tim_datasets
#SBATCH --output=make_all_dataset.log
Rscript $HOME/MS-Thesis/code/networks/make_all_dataset.R
