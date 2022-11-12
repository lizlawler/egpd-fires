#!/bin/bash

module load anaconda
conda activate lawler

Rscript full-model/fire-sims/burns/g2/R/g2_NUTS_sampling.R ${suffix:-NULL} ${params:-NULL}
