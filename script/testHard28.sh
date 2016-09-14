#!/bin/bash

if [ $# -ne 1 ] ; then
echo "Usage: $0 [solver type]"
echo "       solver type 0 is binary search bin packing"
echo "       solver type 1 is recursive number partitioning"
echo "       solver type 2 is belov branch-and-cut-and-price"
exit 1
fi

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


# if no screen, start it
if ! screen -list | grep -q "experiments"; then
  echo "No screen, creating new one with name experiments"
  screen -d -m -S experiments           # Create screen
else
  echo "Screen Exists."
fi

# Binary Search Bin Completion
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem119.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem13.bpa	
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem144.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem14.bpa	
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem175.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem178.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem181.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem195.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem359.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem360.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem40.bpa	
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem419.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem47.bpa	
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem485.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem531.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem561.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem60.bpa	
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem640.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem645.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem709.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem716.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem742.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem766.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem781.bpa

echo "Waiting for a few to finish, be back in 30 seconds"
sleep 30;  # Wait for some to finish

screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem785.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem814.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem832.bpa
screen -S experiments -X screen ./bin/BinPacking -m $1 -f ./dat_hard28/problem900.bpa