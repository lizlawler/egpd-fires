#!/bin/bash

module load anaconda
conda activate lawler

Rscript full-model/fire-sims/burns/g1/R/g1_NUTS_sampling.R ${suffix:-NULL} ${params:-NULL} > full-model/output/g1_$(printf %s ${suffix:-NULL} ${params:-NULL}).txt 2>&1
