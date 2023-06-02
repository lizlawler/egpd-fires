#!/bin/zsh
# model run 

# source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
# conda activate stan

datafile="../../../data/stan_data_${suffix}.json"
basedir="./full-model/fire-sims/${modtype}/${modname}/"
cd ${basedir}
pwd
model="stan/${modname}_${params}"
fitted_files=csv-fits/${modname}_${suffix}_${params}_${delta}

for fit in csv-fits/${modname}_${suffix}_${params}_${delta}*.csv
do
# generate quantities using already fitted parameters
  ./${model} generate_quantities fitted_params=$fit \
                    data file=${datafile} \
                    output file=${fit}_gq.csv &
  echo "Generated quantities for $fit have been calculated"
done
