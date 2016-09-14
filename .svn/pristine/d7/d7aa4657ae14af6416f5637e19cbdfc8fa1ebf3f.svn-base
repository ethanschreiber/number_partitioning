from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})



#
# Function to plot
#

def plot(outFilename,title,minValue,maxValue,allTimes,ALGORITHMS,MARKERS,NArray) :
    with PdfPages(outFilename) as pdf:
        plt.figure(figsize=(5,5))
    
        legend = []

        for i in range(len(ALGORITHMS)) :
            algorithm = ALGORITHMS[i]
            marker = MARKERS[algorithm]
            plt.plot(NArray, allTimes[i],marker,markersize=3,markeredgewidth=1,markerfacecolor='w')
            legend.append(algorithm)
        
        plt.legend(legend, loc='upper left')
        plt.title('Two-Way Partition Results',fontsize=10)
        plt.ylim(minValue,maxValue)
        plt.xlabel('n')
        plt.ylabel('Time (s)')
            
        print("Writing",outFilename)    
        pdf.savefig()
        plt.close()



# ======================
# Read from Command line
# ======================
parser = OptionParser()
parser.add_option("-n", "--minN"  , dest="minN"          ,help="(Required) Min N to process", metavar="int")
parser.add_option("-N", "--maxN"  , dest="maxN"          ,help="(Required) Max N to process", metavar="int")
parser.add_option("-b", "--bits"  , dest="numBits"       ,help="(Required) Number of bits in experiments", metavar="str")
parser.add_option("-x", "--numExp", dest="numExperiments",help="(Required) Num experiments to process", metavar="int")
parser.add_option("-f", "--file"  , dest="filename"      ,help="(Required) Filename to write to", metavar="str")

(options, args) = parser.parse_args()

if (options.minN == None           or options.maxN == None or  
    options.numExperiments == None or options.filename == None) :    
    parser.print_help()
    print ("\n *** Error: Missing required parameters ***\n")
    sys.exit(0)

minN           = int(options.minN)
maxN           = int(options.maxN)
numExperiments = int(options.numExperiments)
numBits        = int(options.numBits)
outFilename       = options.filename


# Make sure parameters are entered
print ("Min N        : " + str(minN))
print ("Max N        : " + str(maxN))
print ("# Experiments: " + str(numExperiments))
print ("Out Filename : " + outFilename)


# ===========
# Run Program
# ===========

if numBits == 16 :
    ALGORITHMS = ["ckk","dp"]
else :
    ALGORITHMS = ["ckk","cga","hs","ss"]

allTimes = []
NArray=[]

MARKERS={"ckk" : "k*--",\
         "cga" : "ks--",\
         "hs"  : "k^--",\
         "ss"  : "ko--"   }
for N in range (minN,maxN+1) : 
   NArray.append(N)

for algorithm in ALGORITHMS :
    times = []
    for N in range (minN,maxN+1) :    
        inputFilename = lib.filename(2,N,algorithm,numExperiments,numBits)
        times.append(lib.processExperiment(inputFilename))

    allTimes.append(times)


maxValue = 150 # max(map(max,allTimes))
minValue = 0 - (maxValue / 100.0)


for times in allTimes :
    for i in range(len(times)) :
        if (times[i] == -1) :
            times[i] = 100000000  # A big number

plot(outFilename, "Optimal Two-Way Partitioning of 48-bit Numbers.",\
     minValue,maxValue,allTimes,ALGORITHMS,MARKERS,NArray)
