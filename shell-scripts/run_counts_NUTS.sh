#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

modtype="counts"
for modname in "zinb" "zinb_er" "zip"
do
for params in "all-reg" "pi-ri"
do
# compile model and link c++ 
object="full-model/fire-sims/${modtype}/${modname}/stan/${modname}_${params}"
cmdstan_model ${object}
for dataset in "climate" "erc" "fwi"
do
sttime=$(date +"%d%b%Y_%H%M")
export modtype modname params dataset sttime
parentjob=$(sbatch --parsable $1 --job-name ${modname}_${params}_${dataset}_${sttime} \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_sampler.sh)
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${modname}_${params}_${dataset}_plots \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_plots.sh
sleep 1
done
done
done