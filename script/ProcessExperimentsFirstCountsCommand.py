import sys
import ProcessExperimentsLibrary

# Paramters for input files
minK = 11
maxK = 12
minN = 40
maxN = 60
NUM_EXPERIMENTS=25

#Parameters for commands
CMD_NUM_EXPERIMENTS=100
MULTIPLIER = 1.5
TIMEOUT = 259200 # 3 days in seconds

lib = ProcessExperimentsLibrary

algorithms = ["ciw","ciw"]

# Data
for N in range (minN,maxN+1) :
    for K in range (minK,maxK+1) :
        count0 = lib.countLines(lib.filename(K,N,algorithms[0],NUM_EXPERIMENTS))
        count1 = lib.countLines(lib.filename(K,N,algorithms[0],NUM_EXPERIMENTS))

        # For Mean
        #t = MULTIPLIER * lib.processExperimentMean(lib.filename(K,N,algorithms[0],NUM_EXPERIMENTS))
        
        # For Max
        t =  MULTIPLIER * lib.processExperimentMax(lib.filename(K,N,algorithms[0],NUM_EXPERIMENTS))

        print "timeout " + str(TIMEOUT) + " ./bin/kPartitioning -f " + \
        "./dat/partitioning/2_48_" + str(N) + "_" + str(CMD_NUM_EXPERIMENTS) + ".np " +\
        " -m 8 -k " + str(K) + " -n " + str(int(t))

        
        
