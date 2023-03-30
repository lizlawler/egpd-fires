#!/bin/bash

#
# modules
#
source /curc/sw/anaconda3/2020.11/etc/profile.d/conda.sh
conda activate /projects/eslawler@colostate.edu/software/anaconda/envs/lawler

#
# run scripts
#
for suffix in "sqrt" "og"
do
for params in "all-reg" "xi-ri"
do
sbatch --job-name g3_$(printf %s $suffix "_" $params) \
--account=csu54_alpine1 \
--chdir=/scratch/alpine/eslawler@colostate.edu/egpd-fires/ \
--output='./full-model/output/%x_%j.txt' --qos=long --nodes=1 --ntasks-per-node=15 \
--time=146:00:00 --mail-type=ALL --mail-user=eslawler@colostate.edu \
--export=suffix=$suffix,params=$params shell-scripts/call_g3_NUTS_sampler.sh
done
done
