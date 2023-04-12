#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

modtype="burns"
modname="g4"
for params in "all-reg" "sigma-ri" "mu-ri"
do
# compile model and link c++ 
object="full-model/fire-sims/${modtype}/${modname}/stan/${modname}_${params}"
cmdstan_model ${object}
for suffix in "sqrt" "og"
do
for delta in 0.81 0.9
do
export modtype modname params suffix delta
parentjob=$(sbatch --parsable $1 --job-name ${modname}_${suffix}_${params}_${delta} \
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