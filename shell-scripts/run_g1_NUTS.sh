#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
module load anaconda
conda activate stan

burn_mod="g1"
for params in "all-reg" "xi-ri" "nu-ri_xi-ri" "kappa-ri_xi-ri" "sigma-ri_xi-ri"
do
# compile model and link c++ 
object="full-model/fire-sims/burns/${burn_mod}/stan/${burn_mod}_${params}"
cmdstan_model ${object}
for suffix in "sqrt" "og"
do
for delta in 0.81 0.9
do
export burn_mod params suffix delta
sbatch --job-name ${burn_mod}_${suffix}_${params}_${delta} \
--output="./full-model/output/%x_%j.txt" \
shell-scripts/call_sampler.sh \
sleep 1
done
done
done