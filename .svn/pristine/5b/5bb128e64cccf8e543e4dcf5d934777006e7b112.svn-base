#!/bin/bash

NUM_FILES=100

#TIME_ALLOWED=172800;  # 2 days in seconds
#TIME_ALLOWED=7200;    # 4 hours in seconds
#TIME_ALLOWED=43200;   # 12 hours in seconds
#TIME_ALLOWED=86400;    # 1 day in seconds
TIME_ALLOWED=172800;   # 2 days in seconds

if [ $# -ne 6 ] ; then
echo ""
echo "   Usage: $0 [min k] [max k] [min n] [max n]  [solver type] [# bits]"
echo "          solver type 0 is binary search bin packing"
echo "          solver type 1 is recursive number partitioning"
echo "          solver type 2 is belov branch-and-cut-and-price"
echo "          solver type 3 is binary search bin packing with Schroeppel and Shamir"
echo ""
exit 1
fi

minK=$1
maxK=$2
minN=$3
maxN=$4
solverType=$5
numBits=$6

DIR=./dat/2_$numBits/

echo ""
echo "Min K      : $minK"
echo "Max K      : $maxK"
echo "Min N      : $minN"
echo "Max N      : $maxN"
echo "Solver Type: $solverType"
echo "Num Bits   : $numBits" 
echo "Output Dir : $DIR"
echo ""

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

# -----------------------------------------
# Create the problems if they do not exist.
# -----------------------------------------
if [ ! -d "$DIR" ]
then
  echo "Creating Directory $PWD"
  ./script/CreateKPartition.sh $numBits $NUM_FILES $minN $maxN
fi

# if no screen, start it
if ! screen -list | grep -q "experiments"; then
  echo "No screen, creating new one with name experiments"
  screen -d -m -S experiments           # Create screen
else
  echo "Screen Exists."
fi

for ((k = $minK; k<=$maxK; k++))      # For each target bin size
do
echo screen -S experiments -X screen ./script/RunKPartition.sh $k $minN $maxN $DIR $solverType $TIME_ALLOWED
     screen -S experiments -X screen ./script/RunKPartition.sh $k $minN $maxN $DIR $solverType $TIME_ALLOWED  
done  

