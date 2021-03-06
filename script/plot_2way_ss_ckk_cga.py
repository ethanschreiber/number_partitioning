from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np

MARKERS={"ckk" : "b-",\
         "cga" : "k-",\
         "ss"  : "r-"
}


ALG_NAMES = {"ckk" : "CKK",\
             "cga" : "CGA",\
             "ss"  : "SS"}
# Plot function
# col 1 is time
# col 3 is isPerfect boolean
# col 4 is memory

def plotExperiments(filename, minBits, maxBits,N, col, ALGORITHMS, minValue, maxValue, isLog): 
    outFile = open(filename,"w")

    allTimes = []
    for i in range(len(ALGORITHMS)) :
        allTimes.append([])
    bitsArray=[]

    for numBits in range (minBits,maxBits+1) :
        bitsArray.append(numBits)

        for i in range(len(ALGORITHMS)) :
            algorithm = ALGORITHMS[i]
            inputFilename = lib.filename(2,N,algorithm,100,numBits)
            print(inputFilename)
            time = lib.processExperiment(inputFilename,100,col)
            allTimes[i].append(time)

    with PdfPages(filename) as pdf:
    
        plt.figure().set_size_inches(7.5,4)
        if (col == 1) :
            plt.title('Average Time to Partition N = ' + str(N)  + ' Integers Two Ways')
            plt.ylabel("Time (s)")
        elif (col == 3) :
            plt.title('Percentage of Perfect Solutions for Instances with n = ' + str(N))
            plt.ylabel("%")
            allTimes[i] = map(lambda(x):x*100,allTimes[i])
        elif (col == 4) :
            plt.title('Average Memory Usage for Instances with n = ' + str(N))
            plt.ylabel("Memory (GB)")
            for i in range(len(allTimes)) :
                allTimes[i] = map(lambda(x):x/1024/1024, allTimes[i])
        else :
            print ("Error, unknow column ", col)
            exit
            
        # See http://matplotlib.org/1.3.1/api/axes_api.html#matplotlib.axes.Axes.plot
        # For styles

        legend = []
        for i in range(len(ALGORITHMS)) :
            algorithm = ALGORITHMS[i]
            marker = MARKERS[algorithm]
            if (col==3) :
                marker = "k--"
            legend.append(ALG_NAMES[algorithm])
            plt.plot(bitsArray,allTimes[i],marker,markersize=5)
            
        if (col != 3):   # No legend for % perfect
            plt.legend(legend, loc='upper left', fontsize=10)
        plt.xlabel("Max # Bits in Input Integers")
        
        if isLog :
            plt.yscale('log')
        plt.grid(True)
        print("Writing",filename)    
        pdf.savefig(bbox_inches='tight')
        plt.close()
    

# ======================
# Read from Command line
# ======================
parser = OptionParser()
parser.add_option("-f", "--file"  , dest="filename"      ,help="(Required) Prefix of filename to write to", metavar="str")

(options, args) = parser.parse_args()

if (options.filename == None) :    
    parser.print_help()
    print ("\n *** Error: Missing required parameters ***\n")
    sys.exit(0)

filename       = options.filename

# Make sure parameters are entered
print ("Out Filename : " + filename)

# ===========
# Run Program
# ===========
maxValue = 400 # max(map(max,allTimes))
minValue = 0 - (maxValue / 100.0)



plotExperiments(filename + "_N_50_runtime.pdf",16,48,50,1,['ckk','cga','ss'],minValue,maxValue, True)
