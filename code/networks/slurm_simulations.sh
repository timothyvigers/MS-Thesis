#!/bin/bash
#SBATCH -n 10
#SBATCH --job-name=tim_simulations
#SBATCH --output=tim_simulations.log
Rscript $HOME/MS-Thesis/code/networks/dic_simulation_simple.R
