#!/bin/bash

source /curc/sw/anaconda3/2020.11/etc/profile.d/conda.sh
conda activate /projects/eslawler@colostate.edu/software/anaconda/envs/lawler

Rscript --vanilla full-model/fire-sims/burns/g3/g3_NUTS_sampling.R \
${suffix:-NULL} ${params:-NULL}
