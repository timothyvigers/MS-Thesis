#!/bin/bash
Rscript step_1/triad_selection_step_1_sv_adj.R
Rscript step_2/triad_selection_step_2_sv_adj.R


Rscript step_1/triad_selection_step_1_sv_unadj.R
Rscript step_2/triad_selection_step_2_sv_unadj.R


Rscript step_1/triad_selection_step_1_long_adj.R