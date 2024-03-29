#!/bin/bash

python ProcessExperiments.py -k 3 -K 7 -n 40 -N 60 -x 100 -f ~/workspace/aaai14/table_3_7.tex 
python ProcessExperiments.py -k 8 -K 12 -n 40 -N 60 -x 100 -f ~/workspace/aaai14/table_8_12.tex 

python ProcessExperimentsSS.py -k 3 -K 7 -n 40 -N 60 -x 100 -f ~/workspace/aaai14/table_SS_3_7.tex 
python ProcessExperimentsSS.py -k 8 -K 12 -n 40 -N 60 -x 100 -f ~/workspace/aaai14/table_SS_8_12.tex 

python plot_hs_ss.py -n 20 -N 60 -x 100 -f ~/workspace/thesis/images/partition
