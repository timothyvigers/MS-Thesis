#!/bin/bash
#SBATCH --time=5
#SBATCH -n 5
#SBATCH --mem=1000
#SBATCH --job-name=tim_slurm_test
#SBATCH --output=tim_slurm_test.log
Rscript candidate_selection/server/step_1/triad_selection_step_1_sv_unadj.R
