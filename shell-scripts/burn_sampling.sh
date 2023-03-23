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

# compile model and link c++ 
# cmdstan_model ${model}
# run model with 3 chains
./${model} sample num_chains=3 num_warmup=1500 num_samples=1500 \
                  adapt delta=${delta} init_buffer=300 term_buffer=200 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv \
                  num_threads=3

# return diagnostics
${diagexe} ${outbase}_*.csv

