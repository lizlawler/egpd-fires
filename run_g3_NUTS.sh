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
sbatch --job-name g3_NUTS_sampling --qos=long --nodes=1 --ntasks-per-node=6 --time=2-00:00:00 --mail-type=ALL --mail-user=eslawler@colostate.edu --export=suffix=$suffix,params=$params call_g3_NUTS_sampler.sh
done
done
