#!/bin/bash
#TIMEOUT=86400  # 1 day in seconds
#TIMEOUT=259200 # 3 days in seconds

minN=20
maxN=70

NUM_PROBLEMS=100


for ((n = $minN; n<=$maxN; n++))                # For each target bin size  
do
    echo "./bin/Partitioning -f ./dat/partitioning/2_48_${n}_$NUM_PROBLEMS.np -m 2 "
    #echo "./bin/Partitioning -f ./dat/partitioning/2_48_${n}_$NUM_PROBLEMS.np -m 3 "
 
done
