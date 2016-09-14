#!/bin/bash

NUM_PROBLEMS=25
B=2   # Base
minE=8 # Exponent
maxE=48

minN=60
maxN=70

NUM_PROBLEMS=25
declare -a METHOD=(1)  # 0=cga, 1=ckk, 2=hs, 3=ss, 4=dp

for e in $(seq $minE 1 $maxE)
do
  for n in $(seq $minN 1 $maxN)
  do
	  for method in ${METHOD[@]}
	  do
  		echo "./bin/Partitioning -f ./dat/partitioning2/2_${e}_${n}_$NUM_PROBLEMS.np -m $method" 
    done
  done
done

