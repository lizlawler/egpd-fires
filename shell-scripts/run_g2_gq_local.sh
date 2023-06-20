#!/bin/zsh
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
trap '' HUP
stanc_exe="/data/accounts/lawler/.conda/envs/stan/bin/cmdstan/bin/stanc"
modtype="burns"
modname="g1"
for params in "all-reg" "xi-ri" "nu-ri_xi-ri" "kappa-ri_xi-ri" "sigma-ri_xi-ri"
do
# compile model and link c++ 
inc_path="full-model/fire-sims/${modtype}/${modname}/stan/"
object="full-model/fire-sims/${modtype}/${modname}/stan/${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for suffix in "og" "sqrt"
do
for delta in 0.81 0.9
do
export modtype modname params suffix delta
nohup ./shell-scripts/burn_gq_mp.sh > full-model/output/${modname}_${suffix}_${params}_${delta}_gq.txt 2>&1 &
sleep 1
done
sleep 1
done
sleep 5000
done
