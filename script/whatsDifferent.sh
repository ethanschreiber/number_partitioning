#!/bin/bash


if [ $# -ne 4 ] ; then
echo ""
echo "   Usage: $0 [k] [min n] [max n] [base dir]"
exit 1
fi

k=$1
minN=$2
maxN=$3
baseDir=$4
  

# ------------------------------------------
# Make sure this is run from BinPacking dir
# ------------------------------------------
PWD=`pwd`
case $PWD in
  *BinPacking ) 
  ;;
  *) echo "Must run from the BinPacking Directory."     
  exit
  ;;    
esac

for ((n = $minN; n<=$maxN; n++))                # For each target bin size  
do


    FILE1=$(printf '%s/%sitems/output_all_%02d_belov_bcp.txt' $baseDir $n $k)
    FILE2=$(printf '%s/%sitems/output_all_%02d_bc.txt' $baseDir $n $k)
     
  ./bin/whatsDifferent $FILE1 $FILE2

done
