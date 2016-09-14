import sys
import os.path

def dataString(time,count=0,NUM_EXPERIMENTS=0) :
 #   count = 0;                   # TEMP TO STOP PRINTING INCOMPLETE COUNTS
 #   NUM_EXPERIMENTS=0;           # TEMP TO STOP PRINTING INCOMPLETE COUNTS
    if (count == 0) :
        return "-"
    # elif (count < NUM_EXPERIMENTS) :
    #     if (time >= 99.9) :
    #         return "%.0f (%d)" % (time,count)
    #     elif (time >= 10) :
    #         return "%.1f (%d)" % (time,count)
    #     elif (time >= 1):
    #         return "%.2f (%d)" % (time,count)
    #     else :
    #         return ("%.3f (%d)" % (time,count)).lstrip('0')
    else :
        if (time >= 99.9) :
            return "%.0f" % time
        elif (time >= 10) :
            return "%.1f" % time
        elif (time >= 1):
            return "%.2f" % time
        else :
            return ("%.2f" % time).lstrip('0')


# For multi-way partitioning
def filename(K, N, algorithm, NUM_EXPERIMENTS,numBits=48) :
    return "/home/ethan/workspace/BinPacking/dat/partitioning/output_" + \
           str(K) + "_2_" + str(numBits) + "_" + str(N) +"_" + \
           str(NUM_EXPERIMENTS) + ".np_" + algorithm + ".txt"

# For two-way partitioning
def filename2(K, N, algorithm, NUM_EXPERIMENTS,numBits=48) :
    return "/home/ethan/workspace/BinPacking/dat/partitioning2/output_" + \
           str(K) + "_2_" + str(numBits) + "_" + str(N) +"_" + \
           str(NUM_EXPERIMENTS) + ".np_" + algorithm + ".txt"

# Return the number of lines in the file or 0
# if the file does not exist
def countLines(filename) :
    try:
        f = open(filename)
        count = 0;
        for line in f:
            count = count + 1           
        return count
    except IOError:
        return 0;

# Read the file and return a list of the times
def processExperimentList(filename, col=1, maxExperiments=1000000) :
    l = []

    
    try:
        f = open(filename)
        
        sum = 0;
        count = 0;
        for line in f:
            lineArr = line.split(" ")
        
            l.append(float(lineArr[col]))
            if (len(l) == maxExperiments) :
                break;
    except IOError:
        print("I/O Error in processExperimentList");
        print "Filename:",filename
    return l;


# TODO: Use processExperimentList
def processExperimentMean(filename) :
    
    try:
        f = open(filename)
        sum = 0;
        count = 0;
        for line in f:
            lineArr = line.split(" ")
            sum += float(lineArr[3])
            count = count + 1
        if (count == 0) :
            return -1
        else :
            return sum / count;
    except IOError:
        return -1;

# TODO: Use processExperimentList
def processExperimentMax(filename) :

    try:
        f = open(filename)
        max = 0;        
        for line in f:
            lineArr = line.split(" ")
            if (float(lineArr[3]) > max) :
                max = float(lineArr[3])

        return max;
    except IOError:
        return -1;

# Only average the first numLines lines of the file
# return 'nan' if file doesn't exist
def processExperiment(filename,numLines=1000000,col=1) :
    
    if (os.path.isfile(filename)) :
        times = processExperimentList(filename,col,numLines)

        if times :                      # if the list is not empty
            return sum(times) / len(times);

    return (float('nan'))


# ==============================================================================
#
# L A T E X   F U N C T I O N S
#
# ==============================================================================

def getAlgorithm(K,i) :
    if (i == 0) :
        return "ciw"
    elif (K <= 7) :
        return "snp"
    elif (K >= 11) :
        return "belov_bcp"
    else :
        return "moffit_rich"

def getName(K,i) :
    if (i == 0) :
        return "CIW"
    elif (K <= 7) :
        return "SNP"
    elif (K >= 11) :
        return "BCP"
    else :
        return "MOF"

def latexTableHeader(minK, maxK, name0=None, name1=None) :
    # Setup tabular
    s = "\\begin{tabular}{r"       # \begin{tabular}{r
    for K in range (minK,maxK+1) :      # |rrr|rrr...
        s += "|rrl"      # print |r for each name
    
    s += "} \\toprule\n"             # }\\ \hline finish line

    # print the first row of header k || 2-way || 3-way ...|
    s += "$k \\rightarrow$ "
    for K in range (minK,maxK+1) :
        if (K == maxK) :
            s += "& \\multicolumn{3}{c}{" + str(K) + "-Way}"
        else :
            s += "& \\multicolumn{3}{c|}{" + str(K) + "-Way}"

        

    s += "\\\\ \n" 

    # print second row of header n | CIW  MOF R | CIW | MOF | R ...
    s += "$n \downarrow$ "
    for K in range (minK,maxK+1) :

        if (name0 == None) :
            s += "& \\multicolumn{1}{c}{" + getName(K,0) + "} "
        else :
            s += "& \\multicolumn{1}{c}{" + name0 + "} "

        if (name1 == None) :
            s += "& \\multicolumn{1}{c}{" + getName(K,1) + "} "
        else :
            s += "& \\multicolumn{1}{c}{" + name1 + "} "

        if (K == maxK) :
            #s += "& \\multicolumn{1}{c}{$\\dfrac{\\text{" + name1 + "}}{\\text{" + name0 + "}}$}"
            s += "& \\multicolumn{1}{c}{R}"
        else :
            #s += "& \\multicolumn{1}{c|}{$\\dfrac{\\text{" + name1 + "}}{\\text{" + name0 + "}}$}"
            s += "& \\multicolumn{1}{c|}{R}"


    s += "\\\\ \\midrule \n"
    return s

def latexTableFooter() :
    s = "\\bottomrule"
    s += "\\end{tabular}"
    return s
