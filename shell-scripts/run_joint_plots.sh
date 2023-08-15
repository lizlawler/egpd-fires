#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate lawler

modtype="joint"
for modname in "all-reg" "kappa-ri" "sigma-ri"
do
for params in "theta-cst_gamma-cst" "theta-cst_gamma-ri" "theta-ri_gamma-cst" "theta-ri_gamma-ri"
do
export modtype modname params
sbatch --job-name ${modname}_${params}_plots \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_plots.sh
sleep 1
done
done
