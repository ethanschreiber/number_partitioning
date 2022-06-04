# number_partitioning

The code is written in C++. It uses some boost libraries and is built with the scons build tool.

# Compiling
1. Install scons: sudo apt-get install scons
2. Install boost: sudo apt-get install libboost-all-dev
From the Binpacking directory, you can build the library using scons:
scons 
The Problem Instances

There are a number of problems for number partitioning defined in dat/partitioning. They are of the format:

2_48_20_100.np

In this case:
2_48 means the input numbers are randomly sampled from the range 0 to 2^48 - 1.
20 means each input problem has 20 numbers.
100 means there are 100 problems in the file.
If you look at a .np file, you will see that each line defines a problem. The first number is the problem number (0 - 99) and the rest of the numbers on the line are the inputs to the problem.

Output Files

Also in dat/partitioning, there are a lot of files beginning with output_. For example:

output_10_2_48_20_100.np_ciw.txt contains the experimental results. In this case:
10 means we partitioned the numbers into 10 subsets.
2_48 means the input numbers are randomly sampled from the range 0 to 2^48-1
20 means each input problem has 20 numbers.
100 means there are 100 problems in the file.
_ciw.txt means we used the cached iterative weakening algorithm. (each suffix should correspond to one of the algorithms from our paper.)

I dont remember what all of the data in the file means, but you could find it somewhere in the source code. The first column is the problem number, the second appears to be the run time, the third is probably the sum of the largest subset. If you can't figure it out from the 
code, I could look deeper to figure out what each of the columns mean.

Command Files and Running

 In the BinPacking directory, there is a script called run_parallel.sh and a  bunch of files starting with command_*
For example, the file command_ciw_11_12.txt defines a set of commands for running cached iterative weakening on various .np files, partitioning into 11 and 12 subsets.

run_parallel is able to solve multiple problem instances at a time, one per cpu core. The algorithm itself is not parallel, but you can run the algorithm on multiple problem instances in parallel.

For example, to run cached iterative weakening on the set of 11-way and 12-way partitioning problems on 10 cores, you can do this:
  

  ./run_parallel.sh command_ciw_11_12.txt 10



If you look through the command files, you will see the problem sets inside of them. You can alternatively run any of the commands in any of the command files directly on the command line. For example:


  ./bin/ciw -f ./dat/partitioning/2_48_20_100.np -k 11

will run cached iterative weakening on 48 bit problems with 20 inputs and partitioning them into 11 subsets.

There are a number of algorithm binaries in the bin directory. If you run any of them with no arguments, there is a help menu.

If there is output for a problem in a dat/partitioning/output_* file, then the problem will be skipped since we already solved it. If you remove dat/partitioning/output_*, you can rerun all of the problem instances from scratch.

Source Files

All of the source files are in src/. Partioning algorithms should all be in src/partition. I have not looked at this code in about six years, so I am not likely to be able to answer in depth  questions, but you can feel free to ask. 
