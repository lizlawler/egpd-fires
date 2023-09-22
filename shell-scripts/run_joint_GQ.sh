#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

stanc_exe="/projects/$USER/software/anaconda/envs/stan/bin/cmdstan/bin/stanc"
modtype="joint"
modname="sigma-ri"
for params in "gamma-ri" "gamma-cst"
do
# compile model and link c++ 
inc_path="full-model/fire-sims/${modtype}/stan/"
object="full-model/fire-sims/${modtype}/stan/${modtype}_${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for iter in 1000 2000
do
for chain in 1 2 3
do
export modtype modname params iter chain
sbatch --job-name ${modname}_${params}_${iter}iter_${chain}_GQ \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_gen_quant.sh
sleep 1
done
done
done