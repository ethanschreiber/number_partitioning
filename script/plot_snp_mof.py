from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np


ALG_NAMES = {"rnp" : "Improved Recursive Number Partitioning",\
             "snp" : "Sequential Number Partitioning",\
             "moffitt" : "Moffitt Partitioning"}
COLOR = {"rnp"     : '1.0',\
         "snp"     : '0.5',\
         "moffitt" : '0.0'}

def plotExperiments(filename, MARKERS, numBits, minN, maxN,  ALGORITHMS): 
    minK=3
    maxK=10
    outFile = open(filename,"w")
    
    NArray=[]
    for N in range (minN,maxN+1) :
        NArray.append(N)
    
    fig, axarr = plt.subplots(4,2, sharex=True)
    fig.set_canvas(plt.gcf().canvas)
    fig.set_size_inches(7.5,9)
    coordsX = 0
    coordsY = 0

    legend = []
    for i in range(len(ALGORITHMS)) :
        legend.append(ALGORITHMS[i])
    for k in range(minK,maxK+1) :
        ax = axarr[coordsX,coordsY]

        allTimes = []
        for i in range(len(ALGORITHMS)) :
            allTimes.append([])
            
        for N in range(minN, maxN+1) :
            for i in range(len(ALGORITHMS)) :
                algorithm = ALGORITHMS[i]
                inputFilename = lib.filename(k,N,algorithm,100,numBits)
                time = lib.processExperiment(inputFilename,100,1)
                allTimes[i].append(time)

                    
        ax.set_title('k = ' + str(k))
                    
        # See http://matplotlib.org/1.3.1/api/axes_api.html#matplotlib.axes.Axes.plot
        # For styles

        for i in range(len(ALGORITHMS)) :
            algorithm = ALGORITHMS[i]
            marker = MARKERS[algorithm]            
            color = COLOR[algorithm]
            #plt.plot(NArray,allTimes[i],marker,markerfacecolor='w',markersize=5)
            #axarr[k-minK].plot(NArray,allTimes[i])
            ax.plot(NArray,allTimes[i],marker,markerfacecolor=color,markersize=7)
        
        if k == minK :
            ax.legend(legend, loc='upper left', fontsize=10)
        if (coordsX == 3) :
            ax.set_xlabel("n")
        
        if (coordsY == 0) :
            ax.set_ylabel("Time (s)")

        ax.set_yscale('log')
        ax.grid(True) 
        #maxValue = max(map(max,allTimes))
        #minValue = 0 - (maxValue / 50.0)
        
        #plt.ylim(minValue,maxValue)

        if coordsY == 1 :
            coordsX = coordsX + 1
            coordsY = 0
        else :
            coordsY = 1

            #plt.suptitle('Average Run Time for 48-Bit Instances')
    print("Writing",filename)    
        #pdf.savefig(bbox_inches='tight')

    fig.savefig(filename, format='pdf',bbox_inches='tight')
    
    

# ======================
# Read from Command line
# ======================
parser = OptionParser()
parser.add_option("-f", "--file"  , dest="filename"      ,help="(Required) Prefix of filename to write to", metavar="str")
parser.add_option("-k", "--k"  , dest="k"      ,help="(Required) Number of subsets.", metavar="int")

(options, args) = parser.parse_args()

if (options.filename == None) :    
    parser.print_help()
    print("\n *** Error: Missing required parameters ***\n")
    sys.exit(0)

filename = options.filename + "_3_10_plot.pdf"

print ("Out Filename : " + filename)


# ===========
# Run Program
# ===========

MARKERS={"rnp" : "ks--",\
         "snp" : "ko--",\
         "moffitt"  : "k^--"\
}

plotExperiments(filename,MARKERS,48,30,45,['snp','moffitt'])
