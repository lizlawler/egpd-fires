#!/bin/bash
# model run 

# change this directory to wherever Stan conda environment lives
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

datafile="../../../data/stan_lists/data_${dataset}.json"
basedir="./models/${modtype}/${modname}/"
cd ${basedir}
model="stan/${modname}_${params}"
outbase="csv_fits/${modname}_${params}_${dataset}_${sttime}"

# run model with 3 chains
./${model} sample num_chains=3 num_warmup=1000 num_samples=1000 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains"