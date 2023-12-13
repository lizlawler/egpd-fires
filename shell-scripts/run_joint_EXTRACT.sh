#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate lawler

modtype="joint"
modname="g1"
params="sigma-ri"
sttime="12Dec2023_1143"
for dataset in "erc" "erc_fwi"
do
export modtype modname params dataset sttime
sbatch --job-name ${modname}_${params}_${sttime}_${dataset}_draws_take2 \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_extraction.sh
sleep 1
sbatch --job-name ${modname}_${params}_${sttime}_${dataset}_scores_take2 \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_plots.sh
sleep 1
done
