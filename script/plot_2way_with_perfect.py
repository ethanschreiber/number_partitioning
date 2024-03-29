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
                    minValue, maxValue, isLog, isGrid, legendLoc, perfectLoc, perfectAlgIdx,mSize): 
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
        
        algorithm = ALGORITHMS[perfectAlgIdx]
        inputFilename = lib.filename(2,N,algorithm,100,numBits)

        perfect = lib.processExperiment(inputFilename,100,3) # 3 is perfect
        perfects.append(perfect * 100)

    with PdfPages(filename) as pdf:
        fig, ax1 = plt.subplots()

        fig.set_size_inches(6.75,2.5)
        ax1.set_title('Average Run Time (left axis) and Percent Perfect Solutions (right axis) \nfor ' + str(numBits) + '-Bit Instances',fontsize=10)        
        ax1.set_ylabel("Time (s)")
            
        # See http://matplotlib.org/1.3.1/api/axes_api.html#matplotlib.axes.Axes.plot
        # For styles

        legend = []
        for i in range(len(ALGORITHMS)) :
            algorithm = ALGORITHMS[i]
            marker = MARKERS[algorithm]
            color = COLOR[algorithm]
        
            legend.append(ALG_NAMES[algorithm])
            
            
            # Change marker size for 32 bit and only show every other marker
            if (numBits == 32) : 
                size = mSize;
                if algorithm == 'cga' :
                    size = 5

                line, = ax1.plot(NArray,allTimes[i],marker,color=color,markersize=size)
                
                line.set_dashes([3,1]) 
            else :
                ax1.plot(NArray,allTimes[i],marker,color=color,markersize=5)
                
            
        ax1.legend(legend, loc=legendLoc, fontsize=10, ncol=4)
        ax1.set_yscale('log')
        

        plt.xlabel("n")
        plt.ylim(minValue,maxValue)
    
        ax2 = ax1.twinx()
        ax2.set_ylabel('% Perfect')
        ax2.set_xlim(20, maxN)
        ax2.set_ylim(-1, 105)

        ax2.plot(NArray,perfects,"--",color="0.0",markersize=5)
        ax2.legend(["perfect %"], loc=perfectLoc, fontsize=10)        
        
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
maxValue48 = 1000000 # max(map(max,allTimes))
minValue48 = 0 - (maxValue48 / 50.0)

maxValue32 = 1000
minValue32 = 0 - (maxValue32 / 50.0)

maxValue16 = 1
minValue16 = 0 - (maxValue16 / 50.0)


MARKERS48={"ckk" : "kD--",\
         "cga" : "k.--",\
         "hs"  : "k^--",\
         "ss"  : "ks--",\
         "dp"  : "ko--"
}


MARKERS32={"ckk" : "kD--",\
           "cga" : "k.--"           
}

MARKERS16={"ckk" : "kD--",\
           "cga" : "k-.",\
           "dp"  : "k-"
}

COLOR = {"ckk": '0.33',\
         "cga": '0.5',\
         "hs" : '0.5',
         "ss" : '0.25',
         "dp" : '0.0'}


COLOR_HS_SS = {"hs" : '1.0',
               "ss" : '0.0'}


plotExperiments(filename + "_32__runtime.pdf",MARKERS48,COLOR,32,20,200,1,['ckk','cga', 'hs', 'ss'],minValue32,maxValue32,True,False,'lower center',"upper right", 0,2)
plotExperiments(filename + "_48__runtime.pdf",MARKERS48,COLOR,48,20,70,1,['ckk','cga', 'hs', 'ss'],minValue48,maxValue48,True,False,'upper left',"lower right", 3,5)
