#!/bin/bash

module load anaconda
conda activate lawler

Rscript full-model/fire-sims/counts/counts_NUTS_sampling.R ${suffix:-NULL} ${params:-NULL} > full-model/output/counts_$(printf %s ${model:-NULL} ${params:-NULL}).txt 2>&1
