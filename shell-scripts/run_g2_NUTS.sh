#!/bin/bash

#
# modules
#
module load anaconda
conda activate lawler

#
# run scripts
#
for suffix in "sqrt" "og"
do
for params in "nu-reg_xi-reg" "nu-reg_xi-ri" "nu-ri_xi-ri"
do
sbatch --job-name g2_$(printf %s $suffix "_" $params) --chdir=/scratch/alpine/eslawler@colostate.edu/egpd-fires/ --output='./full-model/output/%x_%j.txt' --qos=long --nodes=1 --ntasks-per-node=25 --time=72:00:00 --mail-type=ALL --mail-user=eslawler@colostate.edu --export=suffix=$suffix,params=$params shell-scripts/call_g2_NUTS_sampler.sh
done
done
