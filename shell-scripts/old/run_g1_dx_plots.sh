#!/bin/bash

#
# modules
#
module load anaconda
conda activate lawler

#
# run scripts
#
# for suffix in "sqrt" "og"
for suffix in "og"
do
# for params in "all-reg" "xi-ri" "nu-ri_xi-ri" "kappa-ri_xi-ri" "sigma-ri_xi-ri"
for params in "sigma-ri_xi-ri"
do
sbatch --job-name g1_dxplots_$(printf %s $suffix "_" $params) \
--account=csu54_alpine1 \
--chdir=/scratch/alpine/eslawler@colostate.edu/egpd-fires/ \
--output='./full-model/output/%x_%j.txt' --qos=normal --nodes=1 --ntasks-per-node=30 \
--time=1:00:00 --mail-type=ALL --mail-user=eslawler@colostate.edu \
--export=suffix=$suffix,params=$params shell-scripts/call_dx_plots.sh
done
done
