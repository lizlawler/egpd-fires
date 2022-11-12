#!/bin/bash

module load anaconda
conda activate lawler

Rscript full-model/fire-sims/burns/g3/R/g3_NUTS_sampling.R ${suffix:-NULL} ${params:-NULL}
