#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate lawler

modtype="joint"
modname="sigma-ri"
sttime="16Nov2023_1536"
for params in "xi-expit"
do
export modtype modname params sttime
sbatch --job-name ${modname}_${params}_${sttime}_erc_fwi_draws \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_extraction.sh
sleep 1
sbatch --job-name ${modname}_${params}_${sttime}_erc_fwi_scores \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_plots.sh
sleep 1
done
