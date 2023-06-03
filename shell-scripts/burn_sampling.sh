#!/bin/bash
# model run 

source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

datafile="../../../data/stan_data_${suffix}.json"
basedir="./full-model/fire-sims/${modtype}/${modname}/"
# diagexe="/projects/eslawler@colostate.edu/software/anaconda/envs/stan/bin/cmdstan/bin/diagnose"
cd ${basedir}
model="stan/${modname}_${params}"
outbase="csv-fits/${modname}_${suffix}_${params}_${delta}_${sttime}"

# run model with 3 chains
for i in {1..3} 
  do
    ./${model} sample num_warmup=2000 num_samples=2000 thin=2 \
                      adapt delta=${delta} \
                      data file=${datafile} \
                      init=0.01 \
                      output file=${outbase}_${i}.csv &
  done
# # return diagnostics
# ${diagexe} ${outbase}_*.csv
echo "Model has finished running all 3 chains"

