#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

stanc_exe="/projects/$USER/software/anaconda/envs/stan/bin/cmdstan/bin/stanc"
modtype="burns"
modname="g4"
# for params in "all-reg" "xi-ri" "nu-ri_xi-ri" "sigma-ri_xi-ri" "kappa-ri_xi-ri"
for params in "all-reg_cfcns" "xi-ri_cfcns" "nu-ri_xi-ri_cfcns" "sigma-ri_xi-ri_cfcns" "kappa-ri_xi-ri_cfcns"
do
# compile model and link c++ 
inc_path="full-model/fire-sims/${modtype}/${modname}/stan/"
object="full-model/fire-sims/${modtype}/${modname}/stan/${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for suffix in "sqrt" "og"
do
for delta in 0.81 0.9
do
sttime=$(date +"%d%b%Y_%H%M")
export modtype modname params suffix delta sttime
parentjob=$(sbatch --parsable $1 --job-name ${modname}_${suffix}_${params}_${delta}_${sttime} \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_sampler.sh)
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${modname}_${suffix}_${params}_${delta}_plots \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_plots.sh
sleep 1
done
done
done