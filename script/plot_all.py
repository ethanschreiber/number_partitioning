from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np


COLOR = {"cga_mw"   : '1.0',\
         "rnp2009"  : '0.5',\
         "rnp"      : '0.0',\
         "snp"      : '0.5',\
         "moffitt"  : '0.0',\
         "bsbc_IE50": '1.0',\
         "bsbcp"    : '0.0',\
         "ciw"      : '1.0',\
         "ciwlc"    : '0.0'}    
  
def plotExperiments(filename, MARKERS, numBits, minN, maxN, minK, ALGORITHMS,ALG_NAMES): 
    maxK=minK+1
    outFile = open(filename,"w")
    
    NArray=[]
    for N in range (minN,maxN+1) :
        NArray.append(N)
    fig, axarr = plt.subplots(2,1, sharex=False)
    fig.set_canvas(plt.gcf().canvas)
    fig.set_size_inches(7.5,9)

    coordsX = 0

    legend = [] 
    for i in range(len(ALGORITHMS)) :   
        legend.append(ALG_NAMES[ALGORITHMS[i]])

    for k in range(minK,maxK+1) :
        ax = axarr[coordsX]
        ax.grid(True) 
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
                    
        for i in range(len(ALGORITHMS)) :
            algorithm = ALGORITHMS[i]
            marker = MARKERS[algorithm]            
            color = COLOR[algorithm]
            #plt.plot(NArray,allTimes[i],marker,markerfacecolor='w',markersize=5)
            #axarr[k-minK].plot(NArray,allTimes[i])
            ax.plot(NArray,allTimes[i],marker,markerfacecolor=color,markersize=7)
        
        if k == minK :
            ax.legend(legend, loc='lower right', fontsize=10,ncol=2)
        else :
            ax.set_xlabel("n")
        
        ax.set_ylabel("Time (s)")
        ax.set_yscale('log')

        #maxValue = max(map(max,allTimes))
        #minValue = 0 - (maxValue / 50.0)
        
        #plt.ylim(minValue,maxValue)

        coordsX = coordsX + 1

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

# ===========
# Run Program
# ===========


MARKERS = {"cga_mw"   : "k^--",\
           "rnp2009"  : "k^--",\
           "rnp"      : "k^--",\
           "snp"      : "ks--",\
           "moffitt"  : "ks--",\
           "bsbc_IE50": "ko--",\
           "bsbcp"    : "ko--",\
           "ciw"      : "k*--",\
           "ciwlc"    : "k*--"}

ALG_NAMES = {"cga_mw"   : "CGA",\
             "rnp2009"  : "RNP",\
             "rnp"      : "IRNP",\
             "snp"      : "SNPESS",\
             "moffitt"  : "SNPIE",\
             "bsbc_IE50": "BSIBC",\
             "bsbcp"    : "BSBCP",\
             "ciw"      : "CIW",\
             "ciwlc"    : "LCS"}


ALGORITHMS      = ['cga_mw','rnp2009','rnp','snp','moffitt','bsbc_IE50','bsbcp','ciw']

minN=20
maxN=60

for k in range(3,13,2) :
    filename = options.filename + "_" + str(k) + "_" + str(k+1) + "_plot.pdf"
    print ("Out Filename : " + filename)
    plotExperiments(filename,MARKERS,48,minN,maxN,k,ALGORITHMS,ALG_NAMES)


