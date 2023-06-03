#!/bin/zsh
# model run 

trap '' HUP
datafile="../../../data/stan_data_${suffix}.json"
basedir="./full-model/fire-sims/${modtype}/${modname}/"
cd ${basedir}
model="stan/${modname}_${params}"

for fit in csv-fits/${modname}_${suffix}_${params}_${delta}*.csv
do
  gq_file="$(basename -- $fit .csv)"
# generate quantities using already fitted parameters
  nohup ./${model} generate_quantities fitted_params=$fit \
                   data file=${datafile} \
                   output file=csv-fits/${gq_file}_gq.csv &
done