#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

for csv in "g1_og_kappa-ri_xi-ri_0.81_12Apr2023_1322" "g1_og_sigma-ri_xi-ri_0.81_12Apr2023_1323" "g1_og_all-reg_0.81_12Apr2023_1319"
do
export csv
sbatch --job-name ${csv}_plots --output="./full-model/output/%x_%j.txt" \
shell-scripts/call_plots_csv.sh
sleep 1
done