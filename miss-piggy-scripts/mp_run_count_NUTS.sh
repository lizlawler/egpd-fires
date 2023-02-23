#!/bin/zsh

#
# run scripts
#

for model in "zinb" "zinb_er" "zip"
do
for params in "all-reg" "pi-ri"
do
    numjobs='ps -u lawler | grep -e 'g[[:digit:]]' | wc -l'
    while [ $numjobs -gt 6 ]
    do
    echo You have $numbjobs chains running, will fresh in 5 minutes
    sleep 300
    numjobs='ps -u lawler | grep -e 'g[[:digit:]]' | wc -l'
    done

echo Starting model $model with parameters $params
export model=$model params=$params
./shell-scripts/mp_call_count_NUTS_sampler.sh 
done
done