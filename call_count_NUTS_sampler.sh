#!/bin/bash

module load anaconda
conda activate lawler

Rscript full-model/fire-sims/counts/counts_NUTS_sampling.R ${model:-NULL} ${params:-NULL}
