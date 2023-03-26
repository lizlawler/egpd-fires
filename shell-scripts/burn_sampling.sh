#!/bin/bash
# model run 

module load anaconda
conda activate stan

datafile="../../../data/stan_data_${suffix}.json"
basedir="./full-model/fire-sims/burns/${burn_mod}/"
diagexe="/projects/eslawler@colostate.edu/software/anaconda/envs/stan/bin/cmdstan/bin/diagnose"
cd ${basedir}
sttime=$(date +"%d%b%Y_%H%M")
model="stan/${burn_mod}_${params}"
outbase="csv-fits/${burn_mod}_${suffix}_${params}_${delta}_${sttime}"

# run model with 3 chains
./${model} sample num_chains=3 num_warmup=${nwarm} num_samples=2000 thin=2 \
                  adapt delta=${delta} \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv \
                  num_threads=3

# return diagnostics
${diagexe} ${outbase}_*.csv

