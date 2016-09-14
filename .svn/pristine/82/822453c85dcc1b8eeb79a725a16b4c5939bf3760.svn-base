from __future__ import print_function
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

# Read from Command line
parser = OptionParser()
parser.add_option("-k", "--minK"  , dest="minK"          ,help="(Required) Min K to process", metavar="int")
parser.add_option("-K", "--maxK"  , dest="maxK"          ,help="(Required) Max K to process", metavar="int")
parser.add_option("-n", "--minN"  , dest="minN"          ,help="(Required) Min N to process", metavar="int")
parser.add_option("-N", "--maxN"  , dest="maxN"          ,help="(Required) Max N to process", metavar="int")
parser.add_option("-x", "--numExp", dest="numExperiments",help="(Required) Num experiments to process", metavar="int")
parser.add_option("-f", "--file"  , dest="filename"      ,help="(Required) Filename to write to", metavar="str")

(options, args) = parser.parse_args()

minK           = int(options.minK)
maxK           = int(options.maxK)
minN           = int(options.minN)
maxN           = int(options.maxN)
numExperiments = int(options.numExperiments)
filename       = options.filename

# Make sure parameters are entered
if (maxK == None or maxK == None or minN == None or maxN == None or 
    numExperiments == None or filename == None) :    
    parser.print_help()
    print ("\n *** Error: Missing required parameters ***\n")
    sys.exit(0)
else :
    print ("Min K        : " + str(minK))
    print ("Max K        : " + str(maxK))
    print ("Min N        : " + str(minN))
    print ("Max N        : " + str(maxN))
    print ("# Experiments: " + str(numExperiments))
    print ("Out Filename : " + filename)

outFile = open(filename,"w")

print(lib.latexTableHeader(minK,maxK,"CIW","SS"),file=outFile)

# Data
for N in range (minN,maxN+1) :
    print (N,file=outFile,end=''),
    for K in range (minK,maxK+1) :
        algorithm = "ciw"

        count0 = lib.countLines(lib.filename(K,N,algorithm,numExperiments))
        count1 = lib.countLines(lib.filename(K,N,algorithm,numExperiments))
        minCount = min(count0,count1)
        
        t0 = lib.processExperiment(lib.filename(K,N,algorithm,numExperiments),minCount,1)
        t1 = lib.processExperiment(lib.filename(K,N,algorithm,numExperiments),minCount,5)
        ratio = t1 / t0

        if (minCount == 0) : 
            print("& %.1f & - & - " % (t0),file=outFile,end=''),
        else :
            print("& %s & %s & %.2f " % (lib.dataString(t0,count0,numExperiments),lib.dataString(t1,count1,numExperiments),ratio),file=outFile,end=''),
    print("\\\\",file=outFile)
print(lib.latexTableFooter(),file=outFile)

outFile.close()
