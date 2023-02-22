#!/bin/zsh

source /curc/sw/anaconda3/2020.11/etc/profile.d/conda.sh
conda activate /projects/eslawler@colostate.edu/software/anaconda/envs/lawler

Rscript --vanilla full-model/fire-sims/counts/counts_NUTS_sampling.R \
${model:-NULL} ${params:-NULL}
