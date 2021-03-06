from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser
from matplotlib import rc

import matplotlib.pyplot as plt
import numpy as np


COLOR = {"rnp"     : '.5',\
         "rnp2009" : '.5',\
         "snp"     : '0.0',\
         "moffitt" : '1.0',\
         "bsbc_IE50"    : '1.0',\
         "bsbcp"   : '0.0',\
         "ciw"     : '1.0',\
         "ciwlc"   : '0.0'}
   
  
def plotExperiments(filename, MARKERS, numBits, minN, maxN, k, ALGORITHMS,ALG_NAMES): 
    
    pdf = PdfPages(filename) 
    plt.figure(figsize=(8,6))
    plt.subplot(1,1,1)
    plt.yscale('log')
    plt.ylim(.0001,10000)  # TODO: Make theses
    plt.xlim(20,60)        # Paramters
    NArray=[]

    for N in range (minN,maxN+1) :
        NArray.append(N)
        
    allTimes = []
    for i in range(len(ALGORITHMS)) :
        allTimes.append([])
            
    for N in range(minN, maxN+1) :
        for i in range(len(ALGORITHMS)) :
            algorithm = ALGORITHMS[i]
            inputFilename = lib.filename(k,N,algorithm,100,numBits)
            time = lib.processExperiment(inputFilename,100,1)
            allTimes[i].append(time)                          

    rc('text', usetex=True)

    for i in range(len(ALGORITHMS)) :
        algorithm = ALGORITHMS[i]
        marker = MARKERS[algorithm]
        color = COLOR[algorithm]
        #plt.plot(NArray,allTimes[i],marker,markerfacecolor='w',markersize=5)
        #axarr[k-minK].plot(NArray,allTimes[i])
        plt.plot(NArray,allTimes[i],marker,markerfacecolor=color,markersize=7, label=ALG_NAMES[algorithm])
        

    legend = []
    for i in range(len(ALGORITHMS)) :   
        legend.append(ALG_NAMES[ALGORITHMS[i]])
        

    plt.legend( loc='lower right', fontsize=10,ncol=1)
    plt.title('k = ' + str(k) + ' Subsets')
    plt.xlabel("n Input Integers")        
    plt.ylabel("Average Run Time (s)")
    plt.grid(True) 

    pdf.savefig()
    plt.close()
    pdf.close()
    plt.show()

    

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


# MARKERS = {"rnp"      : "ks--",\
#            "snp"      : "ks--",\
#            "moffitt"  : "k+--",\
#            "bsbc_IE50": "ko--",\
#            "bsbcp"    : "ko--",\
#            "ciw"      : "k^--",\
#            "ciwlc"    : "k^--"}

MARKERS = {"rnp"      : "y--",\
           "rnp2009"  : "k--",\
           "snp"      : "g-",\
           "moffitt"  : "r--",\
           "bsbc_IE50": "c-",\
           "bsbcp"    : "m--",\
           "ciw"      : "b-"}

ALG_NAMES = {"rnp"      : "IRNP (2011)",\
             "rnp2009"  : "RNP (2009)" ,\
             "snp"      : r"\textbf{SNP (2013) *}",\
             "moffitt"  : "Moffitt (2013)",\
             "bsbc_IE50": r"\textbf{BSIBC (2013) *}",\
             "bsbcp"    : "BSBCP (2013)",\
             "ciw"      : r"\textbf{CIW (2014) *}",\
             "ciwlc"    : "LCS (2014)"}


ALGORITHMS_MATRIX = [[],[],[],[], # 0 1 2 3
                     ['rnp2009','rnp','snp','ciw'], #4
                     [], [],            # 5 6
                     ['rnp2009','rnp','bsbc_IE50','moffitt','snp','ciw'], #7
                     [], [],            # 8 9
                     ['rnp2009','rnp','bsbc_IE50','bsbcp','moffitt','ciw']] #10
minN=20
maxN=60

for k in range(4,11,3) :
    ALGORITHMS_K = ALGORITHMS_MATRIX[k]
    print("k=" + str(k) + ": ", end="")
    for i in range(len(ALGORITHMS_K)+1) :
        print(str(i) + " ", end="")
        sys.stdout.flush()
        filename = options.filename + "_" + str(k) + "_plot_" + str(i) + ".pdf"
        #print ("Out Filename : " + filename)
        plotExperiments(filename,MARKERS,48,minN,maxN,k,ALGORITHMS_K[0:i],ALG_NAMES)
    print("")

