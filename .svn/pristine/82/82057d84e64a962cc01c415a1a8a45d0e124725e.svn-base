#!/bin/bash
#TIMEOUT=86400  # 1 day in seconds
#TIMEOUT=259200 # 3 days in seconds

declare -a NUM_BITS=(16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60)
declare -a METHOD=(2)  # 0 = CGA, 1=CKK

NUM_PROBLEMS=100
N=50
for numBits in ${NUM_BITS[@]}
do
    for method in ${METHOD[@]}
    do
	echo "./bin/Partitioning -f ./dat/partitioning/2_${numBits}_${N}_$NUM_PROBLEMS.np -m $method" 
    done
done
