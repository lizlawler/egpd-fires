#!/bin/bash
# shell script to run generated_quantities model block

source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

datafile="../../data/stan_data_joint.json"
basedir="./full-model/fire-sims/${modtype}/"
cd ${basedir}
model="stan/${modtype}_${modname}_${params}"

# generate quantities using already fitted parameters
for fit in csv-fits/joint_sigma-ri_${params}_*Sep2023*_${iter}iter_${chain}*
do
  gq_file="../model_comparison/extracted_values/$(basename -- $fit .csv)_redopreds"
  ./${model} generate_quantities fitted_params=$fit \
                            data file=${datafile} \
                          output file=${gq_file}.csv
done

echo "Posterior quantities have been generated"