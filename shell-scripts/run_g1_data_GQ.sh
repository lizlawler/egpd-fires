#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

stanc_exe="/projects/$USER/software/anaconda/envs/stan/bin/cmdstan/bin/stanc"
modtype="burns"
modname="g1"
for params in "sigma-ri_xi-ri"
do
# compile model and link c++ 
inc_path="full-model/fire-sims/${modtype}/${modname}/stan/"
object="full-model/fire-sims/${modtype}/${modname}/stan/${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for dataset in "fwi" "erc_fwi"
do
for qos in "long"
do
for chain in 1 2 3
do
export modtype modname params dataset qos chain
sbatch --job-name ${modname}_${params}_${dataset}_${qos}_${chain}_GQ \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_gen_quant.sh
sleep 1
done
done
done
done