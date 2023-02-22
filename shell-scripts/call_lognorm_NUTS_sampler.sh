#!/bin/bash

source /curc/sw/anaconda3/2020.11/etc/profile.d/conda.sh
conda activate /projects/eslawler@colostate.edu/software/anaconda/envs/lawler

Rscript full-model/fire-sims/burns/lognorm/lognorm_NUTS_sampling.R \
${suffix:-NULL} ${params:-NULL}
