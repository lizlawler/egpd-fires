#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
trap '' HUP
stanc_exe="/data/accounts/lawler/.conda/envs/stan/bin/cmdstan/bin/stanc"
modtype="counts"
modname="zinb_er"
params="pi-ri"
# compile model and link c++ 
inc_path="models/${modtype}/${modname}/stan/"
object="models/${modtype}/${modname}/stan/${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
# for dataset in "erc" "fwi" "erc_fwi"
for dataset in "fwi" "erc_fwi"
do
sttime=$(date +"%d%b%Y_%H%M")
export modtype modname params dataset sttime 
nohup ./shell_scripts/mp_server/submodel_sampling_mp.sh > shell_scripts/console_output/${modname}_${params}_${dataset}_${sttime}_mp.txt 2>&1 &
sleep 1
done