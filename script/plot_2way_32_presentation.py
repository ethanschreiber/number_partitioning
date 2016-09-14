from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np


ALG_NAMES = {"ckk" : "Complete Karmarkar-Karp",\
             "cga" : "Complete Greedy Algorithm",\
             "hs"  : "Horowitz & Sahni",\
             "ss"  : "Schroeppel & Shamir",\
             "dp"  : "Dynamic Programming"}

         

# Plot function
# col 1 is time
# col 3 is isPerfect boolean
# col 4 is memory

def plotExperiments(filename, MARKERS, numBits, minN, maxN, col, ALGORITHMS, 
                    minValue, maxValue, isLog, isGrid, legendLoc): 
    outFile = open(filename,"w")

    allTimes = []
    for i in range(len(ALGORITHMS)) :
        allTimes.append([])
    NArray=[]
    
    perfects = []
    for N in range (minN,maxN+1) :
        NArray.append(N)

        for i in range(len(ALGORITHMS)) :
            algorithm = ALGORITHMS[i]
            inputFilename = lib.filename(2,N,algorithm,100,numBits)
            time = lib.processExperiment(inputFilename,100,col)
            allTimes[i].append(time)
        
        algorithm = ALGORITHMS[0]
        inputFilename = lib.filename(2,N,algorithm,100,numBits)
        perfect = lib.processExperiment(inputFilename,100,3) # 3 is perfect
        perfects.append(perfect * 100)

    with PdfPages(filename) as pdf:
        fig, ax1 = plt.subplots()

        fig.set_size_inches(6.75,5)
        ax1.set_title('Two-Way Partitioning of ' + str(numBits) + '-Bit Instances',fontsize=10)        
        ax1.set_ylabel("Average Run Time (s)")
            
        # See http://matplotlib.org/1.3.1/api/axes_api.html#matplotlib.axes.Axes.plot
        # For styles

        legend = []
        for i in range(len(ALGORITHMS)) :
            algorithm = ALGORITHMS[i]
            marker = MARKERS[algorithm]
        
            legend.append(algorithm)
            ax1.plot(NArray,allTimes[i],marker,markersize=5)
            
        ax1.legend(legend, loc=legendLoc, fontsize=10)        
        ax1.set_yscale('log')
        

        plt.xlabel("n Input Integers")
        plt.ylim(minValue,maxValue)
    
        ax2 = ax1.twinx()
        ax2.set_ylabel('% Perfect')
        ax2.set_xlim(20, 200)
        ax2.set_ylim(-1, 105)

        marker = "g-"
        
        ax2.plot(NArray,perfects,marker)
        ax2.legend(["perfect %"], loc="upper right", fontsize=10)        
        
        algorithm = ALGORITHMS[0]
        marker = MARKERS[algorithm]
        marker="-+"
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
    print("\n *** Error: Missing required parameters ***\n")
    sys.exit(0)

filename       = options.filename
print ("Out Filename : " + filename)

# ===========
# Run Program
# ===========

maxValue32 = 9
minValue32 = 0 - (maxValue32 / 50.0)

MARKERS32={"ckk" : "b-"
}


plotExperiments(filename + "_32__runtime_log.pdf",MARKERS32,32,20,200,1,['ckk'],minValue32,maxValue32,True,False,'lower left')
