#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate lawler

for chain in 1 2 3
do
export chain
sbatch --job-name g1_burn_preds_chain${chain} \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_g1_burn_preds.sh
sleep 1
done