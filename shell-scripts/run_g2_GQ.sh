#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

stanc_exe="/projects/$USER/software/anaconda/envs/stan/bin/cmdstan/bin/stanc"
modtype="burns"
modname="g2"
for params in "sigma-ri_xi-ri" "sigma-ri_xi-ri_nu"
do
# compile model and link c++ 
inc_path="full-model/fire-sims/${modtype}/${modname}/stan/"
object="full-model/fire-sims/${modtype}/${modname}/stan/${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for dataset in "climate" "erc" "fwi" "erc_fwi"
do
for qos in "normal" "long"
do
export modtype modname params dataset qos
sbatch --job-name ${modname}_${params}_${dataset}_${qos}_GQ \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_gen_quant.sh
sleep 1
done
done
done