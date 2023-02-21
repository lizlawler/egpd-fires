#!/bin/zsh

#
# modules
#
source /curc/sw/anaconda3/2020.11/etc/profile.d/conda.sh
conda activate /projects/eslawler@colostate.edu/software/anaconda/envs/lawler

#
# run scripts
#
for model in "zinb" "zinb_er" "zip"
do
for params in "all-reg" "pi-ri"
do
sbatch --job-name count_$(printf %s $model "_" $params) \
--chdir=/scratch/alpine/eslawler@colostate.edu/egpd-fires/ \
--output='./full-model/output/%x_%j.txt' --qos=long --nodes=1 --ntasks-per-node=20 \
--time=72:00:00 --mail-type=ALL --mail-user=eslawler@colostate.edu \
--export=model=$model,params=$params shell-scripts/call_count_NUTS_sampler.sh
done
done
