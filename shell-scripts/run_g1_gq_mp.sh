#!/bin/zsh
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
trap '' HUP
stanc_exe="/data/accounts/lawler/.conda/envs/stan/bin/cmdstan/bin/stanc"
modtype="burns"
modname="g1"
for params in "all-reg_cfcns" "xi-ri_cfcns" "nu-ri_xi-ri_cfcns" "kappa-ri_xi-ri_cfcns" "sigma-ri_xi-ri_cfcns"
do
# compile model and link c++ 
inc_path="full-model/fire-sims/${modtype}/${modname}/stan/"
object="full-model/fire-sims/${modtype}/${modname}/stan/${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for suffix in "sqrt"
do
<<<<<<< HEAD:shell-scripts/run_g1_gq_mp.sh
for delta in 0.81 0.9
=======
for delta in 0.81
>>>>>>> main:shell-scripts/run_g1_GQ.sh
do
export modtype modname params suffix delta
nohup ./shell-scripts/burn_gq_mp.sh > full-model/output/${modname}_${suffix}_${params}_${delta}_gq.txt 2>&1 &
sleep 1
done
done
done
