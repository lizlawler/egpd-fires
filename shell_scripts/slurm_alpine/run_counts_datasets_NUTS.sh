#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

stanc_exe="/projects/$USER/software/anaconda/envs/stan/bin/cmdstan/bin/stanc"
modtype="counts"
modname="zinb_er"
params="pi-ri"
# compile model and link c++ 
inc_path="models/${modtype}/${modname}/stan/"
object="models/${modtype}/${modname}/stan/${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for dataset in "erc" "fwi" "erc_fwi"
do
sttime=$(date +"%d%b%Y_%H%M")
export modtype modname params dataset sttime 
parentjob=$(sbatch --parsable $1 --job-name ${modname}_${params}_${dataset}_${sttime} \
--output="./shell_scripts/console_output/%x_%j.txt" \
shell_scripts/slurm_alpine/call_sampler.sh)
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${modname}_${params}_${dataset}_plots \
--output="./shell_scripts/console_output/%x_%j.txt" \
shell_scripts/slurm_alpine/call_plots_scores.sh
sleep 1
done