#!/bin/bash

NUM_PROBLEMS=100
B=2   # Base
minE=8 # Exponent
maxE=48

minN=20
maxN=100

for e in $(seq $minE 1 $maxE)
do
  for n in $(seq $minN 1 $maxN)
  do
    ./bin/CreatePartitioningProblems -b $B -e $e -z $NUM_PROBLEMS -n $n -s $n
  done
done

