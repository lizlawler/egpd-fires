#!/bin/bash
# model run 

module load anaconda
conda activate stan

datafile="../../../data/stan_data_${suffix}.json"
basedir="./full-model/fire-sims/${modtype}/${modname}/"
diagexe="/projects/eslawler@colostate.edu/software/anaconda/envs/stan/bin/cmdstan/bin/diagnose"
cd ${basedir}
sttime=$(date +"%d%b%Y_%H%M")
model="stan/${modname}_${params}"
outbase="csv-fits/${modname}_${suffix}_${params}_${delta}_${sttime}"

# run model with 3 chains
./${model} sample num_chains=3 num_warmup=1000 num_samples=2000 thin=2 \
                  adapt delta=${delta} \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv \
                  num_threads=3

# return diagnostics
${diagexe} ${outbase}_*.csv

