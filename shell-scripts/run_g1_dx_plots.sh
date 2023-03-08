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
for params in "all-reg" "xi-ri" "nu-ri_xi-ri" "kappa-ri_xi-ri" "sigma-ri_xi-ri"
do
sbatch --job-name g1_dxplots_$(printf %s $suffix "_" $params) \
--account=csu54_alpine1 \
--chdir=/scratch/alpine/eslawler@colostate.edu/egpd-fires/ \
--output='./full-model/output/%x_%j.txt' --qos=normal --nodes=1 --ntasks-per-node=45 \
--time=1:00:00 --mail-type=ALL --mail-user=eslawler@colostate.edu \
--export TMPDIR=/scratch/alpine/eslawler@colostate.edu/ \
--export TMP=/scratch/alpine/eslawler@colostate.edu/ \
--export=suffix=$suffix,params=$params shell-scripts/call_dx_plots.sh
done
done