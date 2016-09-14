#!/bin/bash

NUM_PROBLEMS=100
N=50   # Base

minE=16
maxE=60

for e in $(seq $minE 1 $maxE)
do
  ./bin/CreatePartitioningProblems -b 2 -e $e -z $NUM_PROBLEMS -n $N -s $N
done

#  ./bin/CreatePartitioningProblems -b 2 -e 48 -z 100 -n 40 -s 140

