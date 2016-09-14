# This is for 2-way partitioning
from __future__ import print_function
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser
import subprocess

    
# Read from Command line
parser = OptionParser()
parser.add_option("-b", "--minB"  , dest="minB"          ,help="(Required) Min B to process", metavar="int")
parser.add_option("-B", "--maxB"  , dest="maxB"          ,help="(Required) Max B to process", metavar="int")
parser.add_option("-n", "--N"      , dest="N"            ,help="(Required) Min N to process", metavar="int")
parser.add_option("-x", "--numExp", dest="numExperiments",help="(Required) Num experiments to process", metavar="int")


(options, args) = parser.parse_args()

if (options.minB == None or options.maxB == None or 
    options.N == None or options.numExperiments == None ) :    
    parser.print_help()
    print ("\n *** Error: Missing required parameters ***\n")
    sys.exit(0)

minB           = int(options.minB)
maxB           = int(options.maxB)
N              = int(options.N)
numExperiments = int(options.numExperiments)


# Make sure parameters are entered
print ("Min B        : " + str(minB))
print ("Max B        : " + str(maxB))
print ("N            : " + str(N))
print ("# Experiments: " + str(numExperiments))


ALGORITHMS = ["ckk","ss","cga"]

# Idx from /bin/Partitioning
ALG_IDX   = {"ss"   : "3",\
             "ckk"  : "1",\
             "cga"  : "0"
}
 
HOME = "/home/ethan/workspace/BinPacking/"
FILE_DIR = HOME + "dat/partitioning2/"
minIdx = 0
for B in range (minB,maxB+1) :                                        # for each B value
    bestTime = 3600
    minTimeout = 1
    
    if minIdx > 0 :
        tmp = ALGORITHMS[0]
        ALGORITHMS[0] = ALGORITHMS[minIdx]
        ALGORITHMS[minIdx] = tmp

    for i in range(len(ALGORITHMS)) :

        alg = ALGORITHMS[i]
        
        INPUT_FILENAME = "2_%d_%d_%d.np" %(B, N, numExperiments)
        command        = HOME + "bin/Partitioning"
        inputFilename  = FILE_DIR + INPUT_FILENAME
        outputFilename = FILE_DIR + "output_2_" + INPUT_FILENAME + "_" + alg + ".txt"

        subprocess.call(["timeout",str(max(bestTime,minTimeout)),command,"-m",ALG_IDX[alg],"-f",inputFilename])

        count = lib.countLines(outputFilename)
        if (count == numExperiments) :
            t =( lib.processExperiment(outputFilename,count)+.1) * numExperiments

            if (t < bestTime) :                
                bestTime = ceil(t)
                minIdx = i
                

