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
for params in "xi-expit_v2"
do
# compile model and link c++ 
inc_path="full-model/fire-sims/${modtype}/stan/"
object="full-model/fire-sims/${modtype}/stan/${modtype}_${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
sttime=$(date +"%d%b%Y_%H%M")
export modtype modname params sttime
parentjob=$(sbatch --parsable $1 --job-name ${modname}_${params}_${sttime}_erc_fwi \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_sampler.sh)
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${modname}_${params}_erc_fwi_plots \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_plots.sh
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${modname}_${params}_erc_fwi_draws \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_extraction.sh
done
