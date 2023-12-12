#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

stanc_exe="/projects/$USER/software/anaconda/envs/stan/bin/cmdstan/bin/stanc"
modtype="joint"
modname="g1"
params="sigma-ri"
# compile model and link c++ 
inc_path="full-model/fire-sims/${modtype}/stan/"
object="full-model/fire-sims/${modtype}/stan/${modtype}_${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for dataset in "erc" "erc_fwi"
do
sttime=$(date +"%d%b%Y_%H%M")
export modtype modname params dataset sttime
parentjob=$(sbatch --parsable $1 --job-name ${modtype}_${modname}_${params}_${sttime}_${dataset} \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_sampler.sh)
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${modtype}_${modname}_${params}_${dataset}_plots \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_plots.sh
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${modtype}_${modname}_${params}_${dataset}_draws \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_joint_extraction.sh
done
