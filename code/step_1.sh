#!/bin/bash
# Run scripts with time, low priority, no hangup.
nohup bash -c 'time nice -n19 Rscript code/triad_selection_step_1_long_adj.R' &
nohup bash -c 'time nice -n19 Rscript code/triad_selection_step_1_sv_adj.R' &
nohup bash -c 'time nice -n19 Rscript code/triad_selection_step_1_sv_unadj.R' &