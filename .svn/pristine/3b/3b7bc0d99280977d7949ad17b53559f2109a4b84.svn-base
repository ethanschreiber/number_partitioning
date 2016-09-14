#!/bin/bash

TIME_LIMIT=1000

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

# Binary Search Bin Completion

function run { 
# ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/NU_3_0100_25_1.np 
# ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/NU_3_0100_25_3.np 
# ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/NU_3_0100_25_7.np 
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_1_0050_25_3.np -l 110   -u 111    -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_2_0050_25_3.np -l 1092  -u 1105   -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_2_0100_25_2.np -l 1941  -u 1942   -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0010_05_6.np -l 11385 -u 11575  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0050_10_1.np -l 26233 -u 26234  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0050_10_3.np -l 27764 -u 27765  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0050_10_5.np -l 25296 -u 25297  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0050_10_8.np -l 32266 -u 32267  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0050_25_1.np -l 9517  -u 9688   -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0100_25_0.np -l 21169 -u 21172  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0100_25_1.np -l 17197 -u 17199  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0100_25_2.np -l 21572 -u 21575  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0100_25_3.np -l 20842 -u 20844  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0100_25_4.np -l 20568 -u 20571  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0100_25_5.np -l 20695 -u 20697  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0100_25_6.np -l 20021 -u 20023  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0100_25_7.np -l 19272 -u 19274  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0100_25_8.np -l 20598 -u 20600  -m 0 $1
		timeout $TIME_LIMIT ./bin/kPartitioning -f ./dat_scoop/BPA/hard22/U_3_0100_25_9.np -l 19124 -u 19127  -m 0 $1
}


run "-b 50 --lds"
run "-b 50 "
run "-b 250 --lds"
run "-b 250 "
run "-b 2500 --lds"
run "-b 2500 "
