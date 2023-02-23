#!/bin/zsh
output = 
nohup Rscript --vanilla full-model/fire-sims/counts/counts_NUTS_sampling.R ${model:-NULL} ${params:-NULL}


#!/bin/zsh


nohup ping linuxhint.com > ping.out 2>&1 & 
  outbase = 
  nohup Rscript --vanilla full-model/fire-sims/counts/counts_NUTS_sampling.R \
${model:-NULL} ${params:-NULL} > 
  
  
  nohup nice -n -15 Rscript --vanilla full-model/fire-sims/burns/g1/R/g1_NUTS_sampling_all-reg.R > g1_sqrt_all-reg.txt 2>&1 &
  
  
  ps -u lawler | grep -e 'g[[:digit:]]' | wc -l