#!/bin/bash
# model run 

source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

datafile="../../data/stan_lists/data_${modtype}.json"
basedir="./models/${modtype}/"
cd ${basedir}
model="stan/${modtype}_${modname}_${params}"
outbase="csv_fits/${modtype}_${modname}_${params}_${sttime}"

# run model with 3 chains
./${model} sample num_chains=3 num_warmup=1000 num_samples=1000 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains"
