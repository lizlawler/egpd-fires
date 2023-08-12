#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

stanc_exe="/projects/$USER/software/anaconda/envs/stan/bin/cmdstan/bin/stanc"
modtype="joint"
for modname in "kappa-ri" "sigma-ri"
do
for params in "theta-cst_gamma-cst" "theta-cst_gamma-ri" "theta-ri_gamma-cst" "theta-ri_gamma-ri"
do
# compile model and link c++ 
inc_path="full-model/fire-sims/${modtype}/${modname}/"
object="full-model/fire-sims/${modtype}/${modname}/${modtype}_${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
sttime=$(date +"%d%b%Y_%H%M")
export modtype modname params sttime
parentjob=$(sbatch --parsable $1 --job-name ${modname}_${params}_${sttime} \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_sampler.sh)
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${modname}_${params}_plots \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_plots.sh
sleep 1
done
done
