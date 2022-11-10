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
sbatch --job-name g1_$(printf %s $params "_" $suffix) --output='./full-model/output/output_%x.txt'--nodes=1 --ntasks-per-node=4 --time=24:00:00 --mail-type=ALL --mail-user=eslawler@colostate.edu --export=suffix=$suffix,params=$params call_g1_NUTS_sampler.sh
done
done
