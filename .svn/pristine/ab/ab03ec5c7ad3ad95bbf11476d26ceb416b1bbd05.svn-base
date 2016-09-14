from __future__ import print_function
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

    
# Read from Command line
parser = OptionParser()
parser.add_option("-k", "--minK"  , dest="minK"          ,help="(Required) Min K to process", metavar="int")
parser.add_option("-K", "--maxK"  , dest="maxK"          ,help="(Required) Max K to process", metavar="int")
parser.add_option("-n", "--minN"  , dest="minN"          ,help="(Required) Min N to process", metavar="int")
parser.add_option("-N", "--maxN"  , dest="maxN"          ,help="(Required) Max N to process", metavar="int")
parser.add_option("-x", "--numExp", dest="numExperiments",help="(Required) Num experiments to process", metavar="int")
parser.add_option("-f", "--file"  , dest="filename"      ,help="(Required) Filename to write to", metavar="str")

(options, args) = parser.parse_args()

if (options.minK == None or options.maxK == None or 
    options.minN == None or options.maxN == None or 
    options.numExperiments == None or options.filename == None) :    
    parser.print_help()
    print ("\n *** Error: Missing required parameters ***\n")
    sys.exit(0)

minK           = int(options.minK)
maxK           = int(options.maxK)
minN           = int(options.minN)
maxN           = int(options.maxN)
numExperiments = int(options.numExperiments)
filename       = options.filename

#algorithms = ["moffitt","snp","bsbc_IE50","bsbcp","ciw","ciwlc"]
#headers    = ["MOF","SNP","BSIBC","BSBCP","CIW","LCS"] 
algorithms = ["cga_mw","rnp2009","rnp","moffitt","snp","bsbc_IE50","bsbcp","ciw"]
headers    = ["CGA","RNP","IRNP","SNPIE","SNPESS","BSIBC","BSBCP","CIW"] 
algLen = len(algorithms)

# Make sure parameters are entered
print ("Min K        : " + str(minK))
print ("Max K        : " + str(maxK))
print ("Min N        : " + str(minN))
print ("Max N        : " + str(maxN))
print ("# Experiments: " + str(numExperiments))
print ("Out Filename : " + filename)

outFile = open(filename,"w")


# ------------------------------------------------------------------------------
# Header 
# ------------------------------------------------------------------------------
# Setup tabular
headerStr = "\\begin{tabular}{r"       # \begin{tabular}{r
for K in range (minK,maxK+1) :      # |rrr|rrr...
    headerStr += "|"                        # seperator
    for ignore in algorithms :      # print r for each algorithm
        headerStr += "r"

headerStr += "} \\toprule\n"             # }\\ \hline finish line

# print the first row of header k || 2-way || 3-way ...|
headerStr += "$k \\rightarrow$ \n   "
for K in range (minK,maxK+1) :
    if (K == maxK) :
        headerStr += "& \\multicolumn{" + str(algLen) + "}{c}{" + str(K) + "-Way}"
    else :
        headerStr += "& \\multicolumn{" + str(algLen) + "}{c|}{" + str(K) + "-Way}"



headerStr += "\\\\ \n" 

# print second row of header n | CIW  MOF R | CIW | MOF | R ...
headerStr += "$n \downarrow$"
for K in range (minK,maxK+1) :
    headerStr += "\n   "
    for idx in range(0,algLen) :
        if (idx == algLen-1 and K != maxK) :
            headerStr += "& \\multicolumn{1}{c|}{" + headers[idx] + "} "
        else :
            headerStr += "& \\multicolumn{1}{c}{" + headers[idx] + "} "


headerStr += "\\\\\n\\midrule \n"

print (headerStr,file=outFile,end='')



# Data
for N in range (minN,maxN+1) :
    if (N % 2 == 0) :
        print ("\\rowcolor{Gray}  ",file=outFile,end='');
    else :
        print ("                 ",file=outFile,end='');
    print (N,file=outFile,end=''),
    for K in range (minK,maxK+1) :
        #algorithm0 = lib.getAlgorithm(K,0);
        #algorithm1 = lib.getAlgorithm(K,1);

        for idx in range(0,algLen) :
            count = lib.countLines(lib.filename(K,N,algorithms[idx],numExperiments))
            time  = lib.processExperiment(lib.filename(K,N,algorithms[idx],numExperiments),count)
            print(" & %s" % (lib.dataString(time,count,numExperiments)),file=outFile,end=''),

    print("\\\\",file=outFile)
print(lib.latexTableFooter(),file=outFile)

outFile.close()
