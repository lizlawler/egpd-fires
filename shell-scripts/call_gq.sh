#!/bin/zsh
# 
# #SBATCH --partition=amilan
# #SBATCH --account=csu54_alpine1
# #SBATCH --chdir=/scratch/alpine/eslawler@colostate.edu/egpd-fires/
# #SBATCH --qos=normal
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node=30
# #SBATCH --time=3:00:00
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user=eslawler@colostate.edu
# 
# export TMPDIR=/scratch/alpine/$USER/tmp/
# export TMP=${TMPDIR}
# mkdir -p $TMPDIR
# 
# source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
# conda activate stan

./shell-scripts/burn_gq.sh ${modtype} ${modname} ${suffix} ${params} ${delta} ${sttime} &
