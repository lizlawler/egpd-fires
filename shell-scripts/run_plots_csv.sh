#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

for csv in "g1_og_all-reg_0.81_27Mar2023_1528" "g1_og_all-reg_0.9_27Mar2023_1528" "g1_og_kappa-ri_xi-ri_0.81_27Mar2023_1834" "g1_og_kappa-ri_xi-ri_0.9_27Mar2023_1834" \
"g1_og_nu-ri_xi-ri_0.81_27Mar2023_1712" "g1_og_nu-ri_xi-ri_0.9_27Mar2023_1739" "g1_og_sigma-ri_xi-ri_0.81_27Mar2023_1848" "g1_og_sigma-ri_xi-ri_0.9_27Mar2023_1848" \
"g1_og_xi-ri_0.81_27Mar2023_1557" "g1_og_xi-ri_0.9_27Mar2023_1613" "g1_sqrt_all-reg_0.81_27Mar2023_1525" "g1_sqrt_all-reg_0.9_27Mar2023_1523" "g1_sqrt_kappa-ri_xi-ri_0.81_27Mar2023_1834" \
"g1_sqrt_kappa-ri_xi-ri_0.9_27Mar2023_1834" "g1_sqrt_nu-ri_xi-ri_0.81_27Mar2023_1701" "g1_sqrt_nu-ri_xi-ri_0.9_27Mar2023_1701" "g1_sqrt_sigma-ri_xi-ri_0.81_27Mar2023_1834" \
"g1_sqrt_sigma-ri_xi-ri_0.9_27Mar2023_1834" "g1_sqrt_xi-ri_0.81_27Mar2023_1528" "g1_sqrt_xi-ri_0.9_27Mar2023_1528"  
do
export csv
sbatch --job-name ${csv}_plots --output="./full-model/output/%x_%j.txt" \
shell-scripts/call_plots_csv.sh
sleep 1
done