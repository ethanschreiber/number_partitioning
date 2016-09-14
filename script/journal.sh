#!/bin/sh

python plot_all.py 		  -f ~/workspace/journal/images/all
python plot_snp_cga_mof.py
python plot_2way.py               -f ~/workspace/journal/images/partition
python plot_2way_dp_ckk_cga.py    -f ~/workspace/journal/images/partition_bits
python plot_2way_with_perfect.py  -f ~/workspace/journal/images/partition
python plot_ciw_vs_3_12.py        -f ~/workspace/journal/images/ciw_vs_all.pdf


python ProcessExperimentsAllAlgorithms.py  -k 3 -K 3 -n 20 -N 60 -x 100 -f ~/workspace/journal/tex/all_3.tex
python ProcessExperimentsAllAlgorithms.py  -k 4 -K 4 -n 20 -N 60 -x 100 -f ~/workspace/journal/tex/all_4.tex
python ProcessExperimentsAllAlgorithms.py  -k 5 -K 5 -n 20 -N 60 -x 100 -f ~/workspace/journal/tex/all_5.tex
python ProcessExperimentsAllAlgorithms.py  -k 6 -K 6 -n 20 -N 60 -x 100 -f ~/workspace/journal/tex/all_6.tex
python ProcessExperimentsAllAlgorithms.py  -k 7 -K 7 -n 20 -N 60 -x 100 -f ~/workspace/journal/tex/all_7.tex
python ProcessExperimentsAllAlgorithms.py  -k 8 -K 8 -n 20 -N 60 -x 100 -f ~/workspace/journal/tex/all_8.tex
python ProcessExperimentsAllAlgorithms.py  -k 9 -K 9 -n 20 -N 60 -x 100 -f ~/workspace/journal/tex/all_9.tex
python ProcessExperimentsAllAlgorithms.py  -k 10 -K 10 -n 20 -N 60 -x 100 -f ~/workspace/journal/tex/all_10.tex
python ProcessExperimentsAllAlgorithms.py  -k 11 -K 11 -n 20 -N 60 -x 100 -f ~/workspace/journal/tex/all_11.tex
python ProcessExperimentsAllAlgorithms.py  -k 12 -K 12 -n 20 -N 60 -x 100 -f ~/workspace/journal/tex/all_12.tex
