#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate lawler

modtype="joint"
modname="sigma-ri"
for params in "theta-time_gamma-ri"
do
export modtype modname params
sbatch --job-name ${modname}_${params}_erc_fwi_draws \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_extraction.sh
sleep 1
done
