# This is for 2-way partitioning
from __future__ import print_function
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser




def legend(ALGORITHMS,ALG_NAMES) :
    s = "\\\\" +\
        " & " + ALG_NAMES["ckk"] + "&\\multicolumn{20}{l}{- Complete Karmarkar-Karp}" +\
        " & " + ALG_NAMES["cga"] + "&\\multicolumn{20}{l}{- Complete greedy algorithm}\\\\" +\
        " & " + ALG_NAMES["ss"]  + "&\\multicolumn{20}{l}{- Schroeppel and Shamir}" 
    

    return s
def latexTableHeader(minB, maxB) :
    # Setup tabular
    s = "\\setlength{\\tabcolsep}{0pt}\n \\begin{tabular}{lr"  
    
    for B in range (minB,maxB+1) :      #
        s += "r"      # print |r for each name
    
    s += "}\n"

    s += "   &\\multicolumn{" + str(maxB-minB+1) + "}{c}{Precision $b$}\\\\\n   "
    s += "   \\toprule\n"  
    for B in range (minB,maxB+1) :
        s += "&  \\rotatebox{270}{\\hspace{-4mm} " + str(B) + "}" 
    s += "\\\\"
    # print the first row of header k || 2-way || 3-way ...|
    s += "\\midrule"       
 
    return s

def latexTableFooter(minB, maxB) :
    s = "\\midrule\n&"
    
    for B in range (minB,maxB+1) :
        s += "& \\rotatebox{270}{\\hspace{-3mm}" + str(B) + "}" 
    
        
    s += "\\\\\n\\bottomrule\\\\"
    s += "&&\\multicolumn{" + str(maxB-minB-1) + "}{c}{k = Complete Karmarkar-Karp \\hspace{5mm} s = Schroeppel and Shamir}"
    return s
    
# Read from Command line
parser = OptionParser()
parser.add_option("-b", "--minB"  , dest="minB"          ,help="(Required) Min B to process", metavar="int")
parser.add_option("-B", "--maxB"  , dest="maxB"          ,help="(Required) Max B to process", metavar="int")
parser.add_option("-n", "--minN"  , dest="minN"          ,help="(Required) Min N to process", metavar="int")
parser.add_option("-N", "--maxN"  , dest="maxN"          ,help="(Required) Max N to process", metavar="int")
parser.add_option("-x", "--numExp", dest="numExperiments",help="(Required) Num experiments to process", metavar="int")
parser.add_option("-f", "--file"  , dest="filename"      ,help="(Required) Filename to write to", metavar="str")

(options, args) = parser.parse_args()

if (options.minB == None or options.maxB == None or 
    options.minN == None or options.maxN == None or 
    options.numExperiments == None or options.filename == None) :    
    parser.print_help()
    print ("\n *** Error: Missing required parameters ***\n")
    sys.exit(0)

minB           = int(options.minB)
maxB           = int(options.maxB)
minN           = int(options.minN)
maxN           = int(options.maxN)
numExperiments = int(options.numExperiments)
filename       = options.filename


# Make sure parameters are entered
print ("Min B        : " + str(minB))
print ("Max B        : " + str(maxB))
print ("Min N        : " + str(minN))
print ("Max N        : " + str(maxN))
print ("# Experiments: " + str(numExperiments))
print ("Out Filename : " + filename)

outFile = open(filename,"w")

ALGORITHMS = ["ckk","cga","ss"]

ALG_NAMES = {"ckk"         : "\\cellcolor{gray!60} k",\
             "cga"       : "\\cellcolor{gray!45} g",\
             "ss" : "\\cellcolor{gray!100} s"
}



print(latexTableHeader(minB,maxB),file=outFile)

midN = (maxN+minN-10) / 2
# Data
for N in range (minN,maxN+1) :                                            # for each n
        
    if N == midN :
        print(" \\multirow{10}{*}{\\rotatebox{90}{Input Set Cardinality $n$} \\hspace{1mm}}" ,file=outFile,end=''),                           # Print Yes value
    print (" & " + str(N),file=outFile,end=''),                           # Print N value
    for B in range (minB,maxB+1) :                                        # for each B value

        minTime = float("inf")
        minAlg = "ss"
        for i in range(len(ALGORITHMS)) :
            alg = ALGORITHMS[i]
            filename = lib.filename2(2,N,alg,numExperiments,B)

            count = lib.countLines(filename)
            if (N == 54 and B == 41) :
                print(count," - ",filename)
            if (count > 10) :
                t = lib.processExperiment(filename,count)
                if (t < minTime) :
                    minTime = t
                    minAlg = alg
                    if (N == 54 and B == 41) :
                        print(minAlg," ",minTime)

                               
        print(" & " + ALG_NAMES[minAlg],file=outFile,end='')
    print("\\\\",file=outFile)

print(latexTableFooter(minB,maxB),file=outFile)

#print(legend(ALGORITHMS,ALG_NAMES),file=outFile)
print("\\end{tabular}",file=outFile)


outFile.close()
