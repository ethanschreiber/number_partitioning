#!/bin/bash

if [ $# -ne 1 ] ; then
  echo ""
  echo "   Usage: $0 [N]"
  echo "          N - The number of items per experiment"
  echo ""
  exit 1
fi


NUM_ITEMS=$1    # The number of items per experiment
NUM_FILES=25   # The number of files per exponent
BASE=10         # The base of the capacity
#BASE=2

#TIME_ALLOWED=1800;    # 1/2 hour in seconds
TIME_ALLOWED=86400;    # 1 day in seconds

COMMAND_FILE=command_partitioning_problems_$NUM_ITEMS.txt


declare -a eb_arr=(2 4 6 7 8 9 10 15 20 )  # expected number of bins
declare -a e_arr=(4 6 15)

# -----------------------------------------
# Make sure this is run from BinPacking dir
# -----------------------------------------
source assert_binpacking_dir.inc 

# ----------------------------
# Create problems and run file
# ----------------------------
echo ""
echo "- Creating New Problems in ./dat/partitioning."
echo "- Creating $COMMAND_FILE"
echo ""

rm -f $COMMAND_FILE   # Remove command file

for eb in ${eb_arr[@]}
do
  f=$(echo "scale=2; $eb / ($NUM_ITEMS / 2)" |bc -l | sed 's/^\./0./' | sed 's/1.00/1/' | sed 's/0$//'  )
  for e in ${e_arr[@]}  
  do   
			./bin/CreatePartitioningProblems -z $NUM_FILES -n $NUM_ITEMS -b $BASE -e $e -f $f
			echo "timeout $TIME_ALLOWED ./bin/kPartitioning -f ./dat/partitioning/$NUM_ITEMS/${BASE}_$e/$f/ -m 0" >> $COMMAND_FILE     
#			echo "timeout $TIME_ALLOWED ./bin/kPartitioning -f ./dat/partitioning/$NUM_ITEMS/${BASE}_$e/$f/ -m 0 -c -b 1000000" >> $COMMAND_FILE     
#			echo "timeout $TIME_ALLOWED ./bin/kPartitioning -f ./dat/partitioning/$NUM_ITEMS/${BASE}_$e/$f/ -m 1" >> $COMMAND_FILE     
			echo "timeout $TIME_ALLOWED ./bin/kPartitioning -f ./dat/partitioning/$NUM_ITEMS/${BASE}_$e/$f/ -m 2" >> $COMMAND_FILE     
  done
done

echo ""
echo "  Run using ./run_parallel $COMMAND_FILE [num threads]"
echo ""
