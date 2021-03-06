from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np


ALG_NAMES = {"ckk" : "CKK",\
             "cga" : "CGA",\
             "hs"  : "HS",\
             "ss"  : "SS",\
             "dp"  : "DP"}

         

# Plot function
# col 1 is time
# col 3 is isPerfect boolean
# col 4 is memory

def plotExperiments(filename, MARKERS, COLOR, numBits, minN, maxN, col, ALGORITHMS, 
                    minValue, maxValue, isLog, isGrid, legendLoc): 
    outFile = open(filename,"w")

    allTimes = []
    for i in range(len(ALGORITHMS)) :
        allTimes.append([])
    NArray=[]

    for N in range (minN,maxN+1) :
        NArray.append(N)

        for i in range(len(ALGORITHMS)) :
            algorithm = ALGORITHMS[i]
            inputFilename = lib.filename(2,N,algorithm,100,numBits)
            time = lib.processExperiment(inputFilename,100,col)
            allTimes[i].append(time)

    with PdfPages(filename) as pdf:
            
        plt.figure().set_size_inches(6.75,2.5)
        if (col == 1) :
            plt.title('Average Run Time for ' + str(numBits) + '-Bit Instances')
            plt.ylabel("Time (s)")
        elif (col == 3) :
            plt.title('Percentage of Perfect Solutions for ' +str(numBits) + '-Bit Instances')
            plt.ylabel("%")
            allTimes[i] = map(lambda(x):x*100,allTimes[i])
        elif (col == 4) :
            plt.title('Average Memory Usage for ' + str(numBits) + '-Bit Instances')
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
            color = COLOR[algorithm]
            if (col==3) :
                marker = "k-"                        
            legend.append(ALG_NAMES[algorithm])
            plt.plot(NArray,allTimes[i],marker,markerfacecolor=color,markersize=2)
            
        if (col != 3):   # No legend for % perfect
            plt.legend(legend, loc=legendLoc, fontsize=10)
        
        if isLog :
            plt.yscale('log')

        if isGrid :
            plt.grid(True)

        plt.xlabel("n")
        plt.ylim(minValue,maxValue)
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
maxValue48 = 150 # max(map(max,allTimes))
minValue48 = 0 - (maxValue48 / 50.0)

maxValue32 = 9
minValue32 = 0 - (maxValue32 / 50.0)

maxValue16 = .35
minValue16 = 0 - (maxValue16 / 50.0)


MARKERS={"ckk" : "kD--",\
         "cga" : "ks--",\
         "hs"  : "k^--",\
         "ss"  : "ko--",\
         "dp"  : "ko--"
}


MARKERS32={"ckk" : "kD-",\
           "cga" : "ks--"           
}

MARKERS16={"ckk" : "kD--",\
           "cga" : "ks-.",\
           "dp"  : "k-"
}

COLOR = {"ckk": '.33',\
         "cga": '0.5',\
         "hs" : '0.5',
         "ss" : '0.25',
         "dp" : '0.0'}

COLOR_HS_SS = {"hs" : '1.0',
               "ss" : '0.0'}

#plotExperiments(filename + "_48__runtime.pdf",MARKERS  ,COLOR,48,20,80,1,['ckk','cga','hs','ss'],minValue48,maxValue48,True,True,'upper left')
plotExperiments(filename + "_16__runtime.pdf",MARKERS16,COLOR,16,20,100,1,['dp','cga','ckk'],minValue16,maxValue16,False,False,'upper left')
#plotExperiments(filename + "_48__perfect.pdf",MARKERS  ,COLOR,48,20,80,3,['ss']     ,-5,105,False,False,'upper left')
#plotExperiments(filename + "_16__perfect.pdf",MARKERS16,COLOR,16,20,100,3,['cga'] ,-5,105,False,False,'upper left')
plotExperiments(filename + "_48__memory.pdf" ,MARKERS  ,COLOR_HS_SS,48,20,80,4,['hs','ss'],-.25,13,True,True,'upper left')
