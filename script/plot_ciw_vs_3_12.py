from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np


MARKERS={"ciw"     : "ks--",\
         "snp"     : "ko--",\
         "moffitt" : "k^--",\
         "bsbcp"   : "kd--"
}

ALG_NAMES = {"ciw"     : "CIW",\
             "snp"     : "SNPESS",\
             "moffitt" : "SNPIE",\
             "bsbcp"   : "BCP"
}

COLOR = {"ciw"     : '1.0',\
         "snp"     : '0.67',\
         "moffitt" : '0.33',\
         "bsbcp"   : '0.0'
}

COMPARE_ALG = {3 : "snp",\
               4 : "snp",\
               5 : "snp",\
               6 : "snp",\
               7 : "snp",\
               8 : "moffitt",\
               9 : "moffitt",\
               10 : "moffitt",\
               11 : "bsbcp",\
               12 : "bsbcp"
}

def plotExperiments(filename, numBits, minN, maxN): 
    minK=3
    maxK=12
    outFile = open(filename,"w")
    
    NArray=[]
    for N in range (minN,maxN+1) :
        NArray.append(N)
    
    
    fig, axarr = plt.subplots(5,2, sharex=True)
    fig.set_canvas(plt.gcf().canvas)
    fig.set_size_inches(7.5,9)
    coordsX = 0
    coordsY = 0

    for k in range(minK,maxK+1) :
        legend = []
        legend.append(ALG_NAMES['ciw'])
        legend.append(ALG_NAMES[COMPARE_ALG[k]])

        ax = axarr[coordsX,coordsY]
        
        ALGORITHMS = ['ciw',COMPARE_ALG[k]]
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
            ax.plot(NArray,allTimes[i],marker,markerfacecolor=color,markersize=5)
        
            
        ax.legend(legend, loc='upper left', fontsize=10)
        if (coordsX == 4) :        
            ax.set_xlabel("n")
        
        if (coordsY == 0) :
            ax.set_ylabel("Time (s)")

        # Make labels smaller
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.tick_params(axis='both', which='minor', labelsize=8)
        
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


(options, args) = parser.parse_args()

if (options.filename == None) :    
    parser.print_help()
    print("\n *** Error: Missing required parameters ***\n")
    sys.exit(0)

filename = options.filename

print ("Out Filename : " + filename)


# ===========
# Run Program
# ===========



plotExperiments(filename,48,20,60)
