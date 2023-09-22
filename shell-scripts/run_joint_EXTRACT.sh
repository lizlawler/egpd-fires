#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate lawler

modtype="joint"
modname="sigma-ri"
for params in "theta-time_gamma-ri" "theta-time_gamma-cst"
do
for iter in 1000 2000
do
export modtype modname params iter
sbatch --job-name ${modname}_${params}__${iter}iter_draws \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_extraction.sh
sleep 1
done
done
