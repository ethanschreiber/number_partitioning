#!/bin/bash

if [ $# -ne 4 ] ; then
echo "Usage: $0 [input file] [solver type] [min bins] [max bins]"
echo "       solver type 0 is binary search bin packing"
echo "       solver type 1 is recursive number partitioning"
echo "       solver type 2 is belov branch-and-cut-and-price"
exit 1
fi

RADIUS=0
DIAMETER=10 #$((RADIUS * 2))

SUM=$(cat $1 | (read;read; cat) | awk '{s+=$1} END {print s}')

MIN_CAP=$(($SUM/$4))
MAX_CAP=$(($SUM/$3))

echo $MIN_CAP $MAX_CAP
for ((capacity = $MIN_CAP; capacity<=$MAX_CAP; capacity++))      # For each target bin size
do
  ulimit -t 120
  ./BinPackingRange $1 $2 $capacity
done  

        
 
