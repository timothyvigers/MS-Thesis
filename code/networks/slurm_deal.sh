#!/bin/bash
#SBATCH --job-name=tim_deal_pairs
#SBATCH --output=tim_deal_pairs.log
Rscript $HOME/MS-Thesis/code/networks/deal.R
