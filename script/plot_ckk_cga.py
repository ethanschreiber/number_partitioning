from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np

# ======================
# Read from Command line
# ======================
parser = OptionParser()
parser.add_option("-n", "--minN"  , dest="minN"          ,help="(Required) Min N to process", metavar="int")
parser.add_option("-N", "--maxN"  , dest="maxN"          ,help="(Required) Max N to process", metavar="int")
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
filename       = options.filename

# Make sure parameters are entered
print ("Min N        : " + str(minN))
print ("Max N        : " + str(maxN))
print ("# Experiments: " + str(numExperiments))
print ("Out Filename : " + filename)


# ===========
# Run Program
# ===========
outFile = open(filename,"w")

cgaTimes=[]
ckkTimes=[]
NArray=[]
for N in range (minN,maxN+1) :
    
    cgaFilename = lib.filename(2,N,"cga",numExperiments)
    ckkFilename = lib.filename(2,N,"ckk",numExperiments)

    cgaTimes.append(lib.processExperiment(cgaFilename))
    ckkTimes.append(lib.processExperiment(ckkFilename)) 
    NArray.append(N)
    

with PdfPages(filename) as pdf:
    plt.figure(figsize=(5,5))
    plt.plot(NArray, cgaTimes)
    plt.plot(NArray, ckkTimes)
    plt.legend(['CGA', 'CKK'], loc='upper left')
    plt.title('CKK vs CGA')

    print("Writing",filename)    
    pdf.savefig()
    plt.close()
    
