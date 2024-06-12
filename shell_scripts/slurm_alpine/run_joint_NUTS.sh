#!/bin/bash
# shell script to call sbatch
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

stanc_exe="/projects/$USER/software/anaconda/envs/stan/bin/cmdstan/bin/stanc"
modtype="joint"
modname="g1"
params="sigma-ri"
# compile model and link c++ 
inc_path="models/${modtype}/stan/"
object="models/${modtype}/stan/${modtype}_${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
sttime=$(date +"%d%b%Y_%H%M")
export modtype modname params sttime 
parentjob=$(sbatch --parsable $1 --job-name ${modtype}_${modname}_${params}_${sttime} \
--output="./shell_scripts/console_output/%x_%j.txt" \
shell_scripts/slurm_alpine/call_joint_sampler.sh)
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${modtype}_${modname}_${params}_plots \
--output="./shell_scripts/console_output/%x_%j.txt" \
shell_scripts/slurm_alpine/call_joint_plots.sh
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${modtype}_${modname}_${params}_draws \
--output="./shell_scripts/console_output/%x_%j.txt" \
shell_scripts/slurm_alpine/call_joint_extraction.sh
