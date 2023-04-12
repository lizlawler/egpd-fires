#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=csu54_alpine1
#SBATCH --chdir=/scratch/alpine/eslawler@colostate.edu/egpd-fires/
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=96:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eslawler@colostate.edu

export TMPDIR=/scratch/alpine/$USER/tmp/
export TMP=${TMPDIR}
mkdir -p $TMPDIR

source /curc/sw/anaconda3/2022.10/etc/profile.d/conda.sh
conda activate stan

./shell-scripts/burn_sampling.sh ${modtype} ${modname} ${suffix} ${params} ${delta} ${sttime}