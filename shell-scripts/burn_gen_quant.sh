#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

datafile="../../../data/stan_data_${dataset}.json"
basedir="./full-model/fire-sims/${modtype}/${modname}/"
cd ${basedir}
model="stan/${modname}_${params}"

for fit in csv-fits/${modname}_${params}_${dataset}_*Sep2023*_${qos}*
do
  gq_file="$(basename -- $fit .csv)"
# generate quantities using already fitted parameters
  ./${model} generate_quantities fitted_params=$fit \
                            data file=${datafile} \
                          output file=../../model_comparison/extracted_values/gq_${gq_file}.csv &
done