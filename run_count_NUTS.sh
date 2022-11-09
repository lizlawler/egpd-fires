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
sbatch --job-name count_NUTS_sampling --qos=long --nodes=1 --ntasks-per-node=6 --time=2-00:00:00 --mail-type=ALL --mail-user=eslawler@colostate.edu --export=suffix=$suffix,params=$params call_count_NUTS_sampler.sh
done
done
