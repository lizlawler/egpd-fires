#!/bin/bash

#
# modules
#
module load anaconda
conda activate lawler

#
# run scripts
#
for model in "zinb" "zinb_er" "zip"
do
for params in "pi-reg" "pi-ri"
do
sbatch --job-name count_NUTS_sampling --output='./full-model/output/g1_$(printf %s $model _ $params).txt' --nodes=1 --ntasks-per-node=4 --time=24:00:00 --mail-type=ALL --mail-user=eslawler@colostate.edu --export=model=$model,params=$params call_count_NUTS_sampler.sh
done
done
