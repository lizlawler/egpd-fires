#!/usr/bin/zsh
# shell script to kick off sampling
#
# cycle through loop and launch sampling for each combination
#
source /data/accounts/lawler/.zshrc

conda activate stan_new
trap '' HUP
stanc_exe="/data/accounts/lawler/.conda/envs/stan_new/bin/cmdstan/bin/stanc"
modtype="sizes"
modname="g3"
params="xi-ri"
# compile model and link c++ 
inc_path="models/${modtype}/${modname}/stan/"
object="models/${modtype}/${modname}/stan/${modname}_${params}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
dataset="erc_fwi"
sttime=$(date +"%d%b%Y_%H%M")
export modtype modname params dataset sttime
nohup ./shell_scripts/mp_server/submodel_sampling_mp.sh > shell_scripts/console_output/${modname}_${params}_${dataset}_${sttime}_mp.txt 2>&1 &
