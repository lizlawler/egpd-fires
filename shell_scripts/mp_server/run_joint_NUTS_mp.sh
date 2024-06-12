#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
trap '' HUP
stanc_exe="/data/accounts/lawler/.conda/envs/stan/bin/cmdstan/bin/stanc"
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
nohup ./shell_scripts/mp_server/jointmodel_sampling_mp.sh > shell_scripts/console_output/${modtype}_${modname}_${params}_${sttime}_mp.txt 2>&1 &
sleep 1
done
done
