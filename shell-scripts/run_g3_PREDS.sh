#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate lawler

for chain in 1 3
do
for draw in "less" "greater"
do
export chain draw
sbatch --job-name g3_burn_preds_chain${chain}_${draw}than500 \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_g3_burn_preds.sh
sleep 1
done
done