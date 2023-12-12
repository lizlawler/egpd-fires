#!/bin/bash
# model run 
trap '' HUP
datafile="../../data/stan_data_joint_${dataset}.json"
basedir="./full-model/fire-sims/${modtype}/"
cd ${basedir}
model="stan/${modtype}_${modname}_${params}"
outbase="csv-fits/${modtype}_${modname}_${params}_${sttime}_${dataset}"

# run model with 3 chains
for i in {1..5}
  do
    ./${model} sample num_warmup=1000 num_samples=1000 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}_${i}.csv &
  done
