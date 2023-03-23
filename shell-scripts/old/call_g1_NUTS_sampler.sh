#!/bin/bash

--account=csu54_alpine1 \
--chdir=/scratch/alpine/eslawler@colostate.edu/egpd-fires/ \
--qos=long --nodes=1 --ntasks-per-node=16 \
--time=48:00:00 --mail-type=ALL --mail-user=eslawler@colostate.edu \


Rscript --vanilla full-model/fire-sims/burns/g1/g1_NUTS_sampling.R \
${suffix:-NULL} ${params:-NULL}
