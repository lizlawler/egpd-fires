#!/bin/bash
# shell script to call sbatch
#
# execute all model run shell scripts 
#

shell-scripts/miss-piggy/run_g3_NUTS_mp.sh
sleep 1
shell-scripts/miss-piggy/run_g4_NUTS_mp.sh