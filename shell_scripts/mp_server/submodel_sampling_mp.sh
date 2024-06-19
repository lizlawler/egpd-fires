#!/bin/bash
# model run 
trap '' HUP
datafile="../../../data/stan_lists/data_${dataset}.json"
basedir="./models/${modtype}/${modname}/"
cd ${basedir}
model="stan/${modname}_${params}"
outbase="csv_fits/${modname}_${params}_${dataset}_${sttime}"

# run model with 5 chains
for i in {1..4}
  do
    ./${model} sample num_warmup=1000 num_samples=1000 save_warmup=true \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}_${i}.csv \
                  refresh=25 &
  done

# ./${model} sample num_chains=3 num_warmup=1000 num_samples=1000 \
#                   data file=${datafile} \
#                   init=0.01 \
#                   output file=${outbase}.csv \
#                   num_threads=3

echo "Model has finished running all 3 chains"