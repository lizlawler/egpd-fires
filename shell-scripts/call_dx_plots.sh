#!/bin/bash

source /curc/sw/anaconda3/2020.11/etc/profile.d/conda.sh
conda activate /projects/eslawler@colostate.edu/software/anaconda/envs/lawler

Rscript --vanilla full-model/fire-sims/burns/g1/g1_dx_plots.R \
${suffix:-NULL} ${params:-NULL}
