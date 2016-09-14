from __future__ import print_function
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

numExperiments = 100
ALGORITHMS = ["ERROR","ERROR","ERROR", "snp", "snp", "snp", "snp", "snp", "moffitt", "moffitt", "moffitt", "belov_bcp", "belov_bcp"]

ALG_NAMES = {"snp" : "SNP",\
             "mof" : "MOF",\
             "belov_bcp" : "BCP"}

def printTable(filename,minN,maxN, minK, maxK) :
    outFile = open(filename,"w")
    print(lib.latexTableHeader(minK,maxK),file=outFile)

    # Data
    for N in range (minN,maxN+1) :
#        if (N % 2 == 0) :
#            print ("\\rowcolor{Gray}",file=outFile,end='');
        print (N,file=outFile,end=''),
        for K in range (minK,maxK+1) :
            alg1 = "ciw"
            alg2 = ALGORITHMS[K]

            count0 = lib.countLines(lib.filename(K,N,alg1,numExperiments))
            count1 = lib.countLines(lib.filename(K,N,alg2,numExperiments))
            minCount = min(count0,count1)
            t0 = lib.processExperiment(lib.filename(K,N,alg1,numExperiments),minCount)
            t1 = lib.processExperiment(lib.filename(K,N,alg2,numExperiments),minCount)
        
            if ((count0 == 0) or (count1 == 0) or (t0 == 0) or (t1 == 0)) :
                ratio = "-"
            else :
                if (t1 > t0) :
                    ratio = str(int(round(t1 / t0)))
                else :
                    value = round(t0 / t1)
                    if value > 1 :
                        ratio = "1/" + str(int(value))
                    else :
                        ratio = str(int(value))

            print("& %5s & %5s & %5s " % (lib.dataString(t0,count0,numExperiments),lib.dataString(t1,count1,numExperiments),ratio),file=outFile,end=''),

        print("\\\\",file=outFile)
    print(lib.latexTableFooter(),file=outFile)

    outFile.close()

    
# Read from Command line
parser = OptionParser()
parser.add_option("-k", "--minK"  , dest="minK"          ,help="(Required) Min K to process", metavar="int")
parser.add_option("-K", "--maxK"  , dest="maxK"          ,help="(Required) Max K to process", metavar="int")
parser.add_option("-n", "--minN"  , dest="minN"          ,help="(Required) Min N to process", metavar="int")
parser.add_option("-N", "--maxN"  , dest="maxN"          ,help="(Required) Max N to process", metavar="int")
parser.add_option("-f", "--file"  , dest="filename"      ,help="(Required) Filename to write to", metavar="str")

(options, args) = parser.parse_args()

if (options.minK == None or options.maxK == None or 
    options.minN == None or options.maxN == None or 
    options.filename == None) :    
    parser.print_help()
    print ("\n *** Error: Missing required parameters ***\n")
    sys.exit(0)

minK           = int(options.minK)
maxK           = int(options.maxK)
minN           = int(options.minN)
maxN           = int(options.maxN)
filename       = options.filename



# Make sure parameters are entered
print ("Min K        : " + str(minK))
print ("Max K        : " + str(maxK))
print ("Min N        : " + str(minN))
print ("Max N        : " + str(maxN))
print ("Out Filename : " + filename)





printTable(filename,minN,maxN, minK, maxK)
