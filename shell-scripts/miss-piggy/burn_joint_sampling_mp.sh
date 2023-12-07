#!/bin/bash
# model run 
trap '' HUP
datafile="../../data/stan_data_joint_erc_fwi.json"
basedir="./full-model/fire-sims/${modtype}/"
cd ${basedir}
model="stan/${modtype}_${modname}_${params}"
outbase="csv-fits/${modtype}_${modname}_${params}_${sttime}"

# run model with 3 chains
./${model} sample num_chains=3 num_warmup=1000 num_samples=1000 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains"
