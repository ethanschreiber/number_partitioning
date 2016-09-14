#!/bin/bash

NUM_PROBLEMS=10
B=2   # Base
E=48 # Exponent

minN=20
maxN=200

for n in $(seq $minN 1 $maxN)
do
  ./bin/CreatePartitioningProblems -b $B -e $E -z $NUM_PROBLEMS -n $n -s $n
done

#  ./bin/CreatePartitioningProblems -b 2 -e 48 -z 100 -n 40 -s 140

