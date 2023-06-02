#!/bin/bash
# model run 

source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

datafile="../../../data/stan_data_${suffix}.json"
basedir="./full-model/fire-sims/${modtype}/${modname}/"
cd ${basedir}
model="stan/${modname}_${params}"
fitted="csv-fits/${modname}_${suffix}_${params}_${delta}_${sttime}"

# generate quantities using already fitted parameters
./${model} generate_quantities fitted_params=${fitted}_*.csv \
                  data file=${datafile} \
                  output file=${fitted}_genquant.csv

echo "Generated quantities have all been calculated"