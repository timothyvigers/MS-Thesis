#!/bin/bash
#SBATCH -n 20
#SBATCH --job-name=tim_simulations
#SBATCH --output=tim_simulations.log
Rscript $HOME/MS-Thesis/code/networks/dic_simulation.R
