#!/bin/zsh
ncores=4

#
# modules
#
module load anaconda
conda activate lawler

#
# run scripts
#
for suffix in "sqrt" "og"
do
for params in "nu-reg_xi-reg" "nu-reg_xi-ri" "nu-ri_xi-ri"
do
sbatch --job-name g1_NUTS_sampling --output=$suffix-$params.out --partition=amem --time=3-00:00:00 --ntasks=$ncores --mail-type=ALL --mail-user=eslawler@colostate.edu --export=suffix=$suffix,params=$params call_g1_NUTS_sampler.sh
done
done
