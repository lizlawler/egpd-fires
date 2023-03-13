#!/bin/bash

module load anaconda
conda activate lawler

Rscript --vanilla full-model/fire-sims/burns/g1/g1_dx_plots.R \
${suffix:-NULL} ${params:-NULL}
