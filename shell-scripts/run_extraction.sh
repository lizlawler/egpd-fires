#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

for modtype in "count" "burn"
do
export modtype
sbatch --job-name ${modtype}_extract_gq \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_extract_gq.sh
sleep 1
done