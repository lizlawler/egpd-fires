#!/bin/bash
# model run 
trap '' HUP
datafile="../../../data/stan_data_${suffix}.json"
basedir="./full-model/fire-sims/${modtype}/${modname}/"
cd ${basedir}
model="stan/${modname}_${params}"
outbase="csv-fits/${modname}_${suffix}_${params}_${delta}_${sttime}"

# run model with 3 chains
for i in {1..3} 
do
  nohup ./${model} sample num_warmup=2000 num_samples=2000 thin=2 \
                   adapt delta=${delta} \
                   data file=${datafile} init=0.01 \
                   output file=${outbase}_${i}.csv &
done

