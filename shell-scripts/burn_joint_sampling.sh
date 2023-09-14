#!/bin/bash
# model run 

source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

datafile="../../../data/stan_data_joint.json"
basedir="./full-model/fire-sims/${modtype}/stan/"
cd ${basedir}
model="${modtype}_${modname}_${params}"
outbase="csv-fits/${modtype}_${modname}_${params}_${sttime}_${iter}iter"

# run model with 3 chains
if [ ${iter} -eq 2000 ]
then
  ${thinby}=2
else 
  ${thinby}=1
fi

./${model} sample num_chains=3 num_warmup=${iter} num_samples=${iter} thin=${thinby} \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains"