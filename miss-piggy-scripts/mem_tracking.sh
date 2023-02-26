#!/bin/zsh

njobs=`ps -u lawler | grep -e 'g[[:digit:]]' | wc -l`
while [ $njobs -gt 1 ]
do
  ps -eo pid,command,%mem,time | grep 'g1' >> memory.txt
  echo >> memory.txt 
  sleep 120
  numjobs=`ps -u lawler | grep -e 'g[[:digit:]]' | wc -l`
done
