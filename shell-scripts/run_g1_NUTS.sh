#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
module load anaconda
conda activate stan

modtype="burns"
modname="g1"
for params in "all-reg" "xi-ri" "nu-ri_xi-ri" "kappa-ri_xi-ri" "sigma-ri_xi-ri"
do
# compile model and link c++ 
object="full-model/fire-sims/${modtype}/${modname}/stan/${modname}_${params}"
cmdstan_model ${object}
for suffix in "sqrt" "og"
do
for delta in 0.81 0.9
do
export modtype modname params suffix delta
jobid=$(sbatch --parsable --job-name ${modname}_${suffix}_${params}_${delta} \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_sampler.sh) \
sleep 1
sbatch --dependency=afterok:$jobid \
--job-name ${modname}_${suffix}_${params}_${delta}_plots \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_plots.sh \
sleep 1
done
done
done