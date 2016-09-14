#!/bin/bash

ulimit -s unlimited   # Unlimited stack size

if [ $# -ne 6 ] ; then
echo ""
echo "   Usage: $0 [k] [min n] [max n] [base dir] [solver type] [time allowed]"
echo "          solver type 0 is binary search bin packing"
echo "          solver type 1 is recursive number partitioning"
echo "          solver type 2 is belov branch-and-cut-and-price"
echo "          solver type 3 is binary search bin packing with Schroeppel and Shamir"
echo ""

exit 1
fi

k=$1
minN=$2
maxN=$3
baseDir=$4
solverType=$5
TIME_ALLOWED=$6


for ((n = $minN; n<=$maxN; n++))                # For each target bin size  
do
    DIR="$baseDir${n}items/"
    echo "Dir: $DIR"
    
    timeout $TIME_ALLOWED ./bin/kPartitioning -f $DIR -k $k -m $solverType      
done
