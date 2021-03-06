from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np

# Plot function
# col 1 is time
# col 3 is isPerfect boolean
# col 4 is memory

def plotExperiments(filename, minN, maxN, col): 
    outFile = open(filename,"w")

    hsValues=[]
    ssValues=[]
    cgaValues=[]
    ckkValues=[]
    NArray=[]
    for N in range (minN,maxN+1) :
    
        hsFilename  = lib.filename(2,N,"hs",numExperiments)
        ssFilename  = lib.filename(2,N,"ss",numExperiments)

        hsValues.append(lib.processExperiment(hsFilename,100, col))
        ssValues.append(lib.processExperiment(ssFilename,100, col)) 
        NArray.append(N)
    
        if (col == 1) :
            cgaFilename = lib.filename(2,N,"cga",numExperiments)
            ckkFilename = lib.filename(2,N,"ckk",numExperiments)

            print("N: ",N," CGA: ", lib.processExperiment(cgaFilename,100, col))
            cgaTime = lib.processExperiment(cgaFilename,100, col)
            ckkTime = lib.processExperiment(ckkFilename,100, col)
            if (cgaTime != -1) :
                cgaValues.append(cgaTime)
            else :
                cgaValues.append(float('nan'))

            if (ckkTime != -1) :
                ckkValues.append(ckkTime) 
            else :
                ckkValues.append(float('nan'))

    with PdfPages(filename) as pdf:
        
        plt.figure(figsize=(5,5))            
        plt.ylabel("Time (s)")
        if (col == 1) :
            plt.title('Average Run Time')
            plt.ylabel("Time (s)")            
        elif (col == 3) :
            plt.title('Percentage of Perfect Solutions')
            plt.ylabel("%")
            hsValues = map(lambda(x):x*100, hsValues)
        elif (col == 4) :
            plt.title('Average Memory Usage')
            plt.ylabel("Memory (GB)")
            hsValues = map(lambda(x):x/1024/1024, hsValues)
            ssValues = map(lambda(x):x/1024/1024, ssValues)
        else :
            print ("Error, unknow column ", col)
            exit
            
        # See http://matplotlib.org/1.3.1/api/axes_api.html#matplotlib.axes.Axes.plot
        # For styles
        if (col == 3) :
            plt.plot(NArray, hsValues, 'k-')
        else :
            plt.plot(NArray, hsValues, 'k-')
            plt.plot(NArray, ssValues, 'k--')

            if (col == 1) :
                plt.plot(NArray, cgaValues, 'k:')
                plt.plot(NArray, ckkValues, 'k,')
                plt.legend(['Horowitz & Sahni', 'Schroeppel & Shamir', 'Complete Greedy Algorithm', 'Complete Karmarkar-Karp',], loc='upper left', fontsize=10)
            else :
                plt.legend(['Horowitz & Sahni', 'Schroeppel & Shamir'], loc='upper left', fontsize=10)
            
        plt.xlabel("N")

        plt.ylim(min(hsValues)-1,max(hsValues)+1)
        
        print("Writing",filename)    
        pdf.savefig()
        plt.close()
    

# ======================
# Read from Command line
# ======================
parser = OptionParser()
parser.add_option("-n", "--minN"  , dest="minN"          ,help="(Required) Min N to process", metavar="int")
parser.add_option("-N", "--maxN"  , dest="maxN"          ,help="(Required) Max N to process", metavar="int")
parser.add_option("-x", "--numExp", dest="numExperiments",help="(Required) Num experiments to process", metavar="int")
parser.add_option("-f", "--file"  , dest="filename"      ,help="(Required) Prefix of filename to write to", metavar="str")

(options, args) = parser.parse_args()

if (options.minN == None           or options.maxN == None or  
    options.numExperiments == None or options.filename == None) :    
    parser.print_help()
    print ("\n *** Error: Missing required parameters ***\n")
    sys.exit(0)

minN           = int(options.minN)
maxN           = int(options.maxN)
numExperiments = int(options.numExperiments)
filename       = options.filename

# Make sure parameters are entered
print ("Min N        : " + str(minN))
print ("Max N        : " + str(maxN))
print ("# Experiments: " + str(numExperiments))
print ("Out Filename : " + filename)


# ===========
# Run Program
# ===========
plotExperiments(filename + "_runtime.pdf",minN,maxN,1)
plotExperiments(filename + "_perfect.pdf",minN,maxN,3)
plotExperiments(filename + "_memory.pdf",minN,maxN,4)
