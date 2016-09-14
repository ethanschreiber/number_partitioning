from __future__ import print_function
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser



def legend(ALGORITHMS,ALG_NAMES) :
    s = "\\\\" +\
        " & " + ALG_NAMES["ciw"]       + "&\\multicolumn{4}{l}{- Cached Iterative Weakening }"                + \
        " & " + ALG_NAMES["bsbc_IE50"] + "&\\multicolumn{4}{l}{- Binary-Search Improved Bin Completion} \\\\" + \
        " & " + ALG_NAMES["snp"]       + "&\\multicolumn{4}{l}{- Sequential Number Partitioning}"             + \
        " & " + ALG_NAMES["bsbcp"]     + "&\\multicolumn{4}{l}{- Binary-Search Branch-and-Cut-and-Price}\\\\" + \
        " & " + ALG_NAMES["moffitt"]   + "&\\multicolumn{4}{l}{- Moffitt Partitioning}"                       + \
        " & " + ""                     + "&\\multicolumn{4}{l}{  \;\; (previous work)                      }\\\\" 
        

    

    return s
def latexTableHeader(minK, maxK) :
    # Setup tabular
    s = "\\begin{tabular}{r"       # \begin{tabular}{r
    for K in range (minK,maxK+1) :      #
        s += "l"      # print |r for each name
    
    s += "} \\toprule\n $k \\rightarrow$\\\\"             # }\\ \hline finish line

    # print the first row of header k || 2-way || 3-way ...|
    s += "$n \downarrow$  "
    for K in range (minK,maxK+1) :
        s += "& " + str(K) + "-Way"

    s += "\\\\ \\midrule"       
 
    return s


    
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


# Make sure parameters are entered
print ("Min K        : " + str(minK))
print ("Max K        : " + str(maxK))
print ("Min N        : " + str(minN))
print ("Max N        : " + str(maxN))
print ("# Experiments: " + str(numExperiments))
print ("Out Filename : " + filename)

outFile = open(filename,"w")

ALGORITHMS = ["ciw","moffitt","bsbc_IE50","bsbcp","rnp","snp"]

ALG_NAMES = {"ciw"         : "\\cellcolor{gray!85} CIW",\
             "snp"         : "\\cellcolor{gray!70} SNP",\
             "bsbc_IE50"   : "\\cellcolor{gray!55} IBC",\
             "moffitt"     : "\\cellcolor{gray!40} MOF",\
             "bsbcp"       : "\\cellcolor{gray!25} BCP"}             




print(latexTableHeader(minK,maxK),file=outFile)

# Data
for N in range (minN,maxN+1) :                                            # for each n
        
    print ("   " + str(N),file=outFile,end=''),                                        # Print N value
    for K in range (minK,maxK+1) :                                        # for each K value

        minTime = float("inf")
        minAlg = "ciw"
        for i in range(len(ALGORITHMS)) :
            alg = ALGORITHMS[i]
            count = lib.countLines(lib.filename(K,N,alg,numExperiments))
            if (count > 0) :
                t = lib.processExperiment(lib.filename(K,N,alg,numExperiments),count)

                if (t < minTime) :
                    minTime = t
                    minAlg = alg
                                
        print(" & " + ALG_NAMES[minAlg],file=outFile,end='')
    print("\\\\",file=outFile)

print("\\bottomrule",file=outFile)

print(legend(ALGORITHMS,ALG_NAMES),file=outFile)
print("\\end{tabular}",file=outFile)


outFile.close()
