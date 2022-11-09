#!/bin/zsh

module load anaconda
conda activate lawler

Rscript ./full-model/fire-sims/burns/g1/R/g1_NUTS_sampling.R ${suffix:-NULL} ${params:-NULL}