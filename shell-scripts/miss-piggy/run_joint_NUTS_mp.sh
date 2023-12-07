#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
trap '' HUP
stanc_exe="/data/accounts/lawler/.conda/envs/stan/bin/cmdstan/bin/stanc"
modtype="joint"
for modname in "g3" "g4"
do
for params in "nu" "sigma"
do
# compile model and link c++ 
inc_path="full-model/fire-sims/${modtype}/${modname}/stan/"
object="full-model/fire-sims/${modtype}/${modname}/stan/${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
sttime=$(date +"%d%b%Y_%H%M")
export modtype modname params sttime
nohup ./shell-scripts/miss-piggy/burn_joint_sampling_mp.sh > full-model/output/${modname}_${dataset}_${params}_${sttime}_mp.txt 2>&1 &
sleep 1
done
done
