#!/bin/bash
#TIMEOUT=86400  # 1 day in seconds
#TIMEOUT=259200 # 3 days in seconds


minN=20
maxN=30
minK=10
maxK=10
#declare -a K=(3 4 5 6 7)              # For SNP
#declare -a K=(11 12)                   # For BSBCP
#declare -a K=(8 9 10)                 # For MOF

declare -a METHOD=(rnp2009)

NUM_PROBLEMS=100
BITS=48

for ((n = $minN; n<=$maxN; n++))      
do
    for ((k = $minK; k<=$maxK; k++))  
    do
  	  for method in ${METHOD[@]}
	    do
	      echo "./bin/$method -f ./dat/partitioning/2_${BITS}_${n}_${NUM_PROBLEMS}.np -k $k" 
	    done
    done
done
