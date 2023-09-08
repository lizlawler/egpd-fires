#!/bin/bash
# shell script to call sbatch
#
# execute all model run shell scripts 
#
source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

shell-scripts/run_joint_NUTS.sh
sleep 1
shell-scripts/run_g1_g2_NUTS.sh
sleep 1
shell-scripts/run_counts_NUTS.sh
sleep 1
shell-scripts/run_lognorm_NUTS.sh
sleep 1
shell-scripts/run_g3_NUTS.sh
sleep 1
shell-scripts/run_g4_NUTS.sh