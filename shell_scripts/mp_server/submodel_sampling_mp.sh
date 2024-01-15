#!/bin/bash
# model run 
trap '' HUP
datafile="../../../data/stan_lists/data_${dataset}.json"
basedir="./models/${modtype}/${modname}/"
cd ${basedir}
model="stan/${modname}_${params}"
outbase="csv_fits/${modname}_${params}_${dataset}_${sttime}"

# run model with 3 chains
for i in {1..4}
  do
    ./${model} sample num_warmup=1000 num_samples=1000 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}_${i}.csv &
  done
