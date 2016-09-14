#!/bin/bash

if [ $# -ne 1 ] ; then
  echo ""
  echo "   Usage: $0 [N]"
  echo "          N - The number of items per experiment"
  echo ""
  exit 1
fi


NUM_ITEMS=$1    # The number of items per experiment
NUM_FILES=100   # The number of files per exponent
BASE=10    # The base of the capacity

#TIME_ALLOWED=1800;    # 1/2 hour in seconds
#TIME_ALLOWED=86400;    # 1 day in seconds
TIME_ALLOWED=21600;    # 6 hours in seconds

COMMAND_FILE=command_packing_problems_$NUM_ITEMS.txt

declare -a eb_arr=(2 4 6 7 8 9 10 15 20 21 22 23 24 25)  # expected number of bins for packing, exact for partitioning
declare -a e_arr=(4 6 15)                            # exponent. i.e. 10^4 10^6, 10^15

# -----------------------------------------
# Make sure this is run from BinPacking dir
# -----------------------------------------
source assert_binpacking_dir.inc 

# ----------------------------
# Create problems and run file
# ----------------------------
echo ""
echo "- Creating New Problems in ./dat/packing."
echo "- Creating $COMMAND_FILE"
echo ""

rm -f $COMMAND_FILE   # Remove command file

for eb in ${eb_arr[@]}
do
  f=$(echo "scale=2; $eb / ($NUM_ITEMS / 2)" |bc -l | sed 's/^\./0./' | sed 's/1.00/1/' | sed 's/0$//'  )
  for e in ${e_arr[@]}  
  do
			echo ""
			./bin/CreatePackingProblems -z $NUM_FILES -n $NUM_ITEMS -b $BASE -e $e -f $f
#			echo "timeout $TIME_ALLOWED ./bin/BinPacking -f ./dat/packing/$NUM_ITEMS/${BASE}_$e/$f/" >> $COMMAND_FILE -m 0    
			echo "timeout $TIME_ALLOWED ./bin/BinPacking -f ./dat/packing/$NUM_ITEMS/${BASE}_$e/$f/" >> $COMMAND_FILE -m 0 -c -b 1000000 
			echo "timeout $TIME_ALLOWED ./bin/BinPacking -f ./dat/packing/$NUM_ITEMS/${BASE}_$e/$f/" >> $COMMAND_FILE -m 0 -b 50
#			echo "timeout $TIME_ALLOWED ./bin/BinPacking -f ./dat/packing/$NUM_ITEMS/${BASE}_$e/$f/" >> $COMMAND_FILE -m 2            
  done
done

echo ""
echo "  Run using ./run_all $COMMAND_FILE"
echo ""
