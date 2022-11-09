#!/bin/bash

module load anaconda
conda activate lawler

Rscript full-model/fire-sims/burns/g2/R/g2_NUTS_sampling.R ${suffix:-NULL} ${params:-NULL} > full-model/output/g2_$(printf %s ${suffix:-NULL} ${params:-NULL}).txt 2>&1
