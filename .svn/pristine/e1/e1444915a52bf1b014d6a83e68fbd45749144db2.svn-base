#!/bin/bash

NUM_FILES=100
MIN=0

MAX=10000
DIR=dat_10000

MAX=281474976710655


if [ $# -ne 4 ] ; then
echo ""
echo "   Usage: $0 [# bits] [# files] [min n] [max n] "
echo "          - All values in the range [0,2^{# bits}]"
echo "          - Generate [# files] for each n."
echo "          - Generate for n from [min n] to [max n]"
echo ""
exit 1
fi

numBits=$1
numFiles=$2
minN=$3
maxN=$4
DIR=./dat/2_$numBits/

echo ""
echo "Num Bits : $numBits"
echo "Num Files: $numFiles"
echo "Min N    : $minN"
echo "Max N    : $maxN"
echo ""

# ------------------------------------------
# Make sure this is run from BinPacking dir
# ------------------------------------------
PWD=`pwd`
case $PWD in
  *BinPacking ) 
  ;;
  *) echo "Must run from the BinPacking Directory: ./script/CreateKPartition.sh"
  exit
  ;;    
esac

# ----------------------------------------------------
# Create the dir if it does not exist. Exit otherwise
# ----------------------------------------------------
if [ ! -d "$DIR" ]
then
  echo "Creating Directory: $DIR"
  mkdir -p $DIR
else
  echo "$DIR exists, exiting!"
  exit;
fi  



for ((n = $minN; n<=$maxN; n++))      # For each target bin size
do
  SEED=$n     # For now, use n as the seed, why not?
    
  ./bin/CreatePartitioningProblems $NUM_FILES $n $SEED $numBits $DIR
done

  
