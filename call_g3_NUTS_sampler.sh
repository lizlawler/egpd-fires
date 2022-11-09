#!/bin/bash

module load anaconda
conda activate lawler

Rscript full-model/fire-sims/burns/g3/R/g3_NUTS_sampling.R ${suffix:-NULL} ${params:-NULL} > full-model/output/g3_$(printf %s ${suffix:-NULL} ${params:-NULL}).txt 2>&1
