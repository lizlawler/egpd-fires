#!/bin/bash
# model run 
trap '' HUP
datafile="../../../data/stan_data_${dataset}.json"
basedir="./full-model/fire-sims/${modtype}/${modname}/"
cd ${basedir}
model="stan/${modname}_${params}"
outbase="csv-fits/${modname}_${dataset}_${params}_${sttime}"

# run model with 3 chains
./${model} sample num_chains=3 num_warmup=1000 num_samples=1000 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains"
echo "Now running script to create trace plots and pull scores..."
export modtype modname params dataset sttime
nohup Rscript --vanilla full-model/fire-sims/dx_plots_scores_mp.R > \
full-model/output/${modname}_${dataset}_${params}_plots_mp.txt 2>&1 &
