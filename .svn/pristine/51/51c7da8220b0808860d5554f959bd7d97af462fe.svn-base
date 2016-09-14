import sys

FIELD_WIDTH = 8
minK = 3
maxK = 10
minN = 40
maxN = 60
NUM_EXPERIMENTS=100

def dataString(time,count) :
    if (count < NUM_EXPERIMENTS) :
        return "%d (%d)" % (time,count)
    else :
        return "%d" % time

def filename(K,N,algorithm) :
    return "/home/ethan/workspace/BinPacking/dat/partitioning/output_" + str(K) + "_2_48_" + str(N) +"_" + str(NUM_EXPERIMENTS) + ".np_" + algorithm + ".txt"


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


algorithms = ["ciw","ciw"]

# Setup tabular
print("\\begin{tabular}{|r"),       # \begin{tabular}{|r
for K in range (minK,maxK+1) :      # ||r|r||r|r...
    sys.stdout.write("||r|r")           # Start with double || so print first |
    
print("|} \\hline")             # |}\\ \hline finish line

# print the first row of header k || 2-way || 3-way ...|
print("\\multicolumn{1}{|c||}{$k$} "),
for K in range (minK,maxK+1) :
    if (K == maxK) :
        print("& \\multicolumn{2}{|c|}{" + str(K) + "-Way}"),
    else :
        print("& \\multicolumn{2}{|c||}{" + str(K) + "-Way}"),
print("\\\\ \\hline")        

# print second row of header n || IP | MOF | R || IP | MOF | R ...|
print("\\multicolumn{1}{|c||}{$n$} ")
for K in range (minK,maxK+1) :

    print("& \\multicolumn{1}{|c}{Mean} "),
    if (K == maxK) :
        print("& \\multicolumn{1}{|c|}{Max} "),
    else :
        print("& \\multicolumn{1}{|c||}{Max} "),


print("\\\\ \\hline")        

# Data
for N in range (minN,maxN+1) :
    print (N),
    for K in range (minK,maxK+1) :
        count0 = countLines(filename(K,N,algorithms[0]))
        count1 = countLines(filename(K,N,algorithms[0]))
        
        t0 = processExperimentMean(filename(K,N,algorithms[0]))
        t1 = processExperimentMax(filename(K,N,algorithms[0]))
        print("& %s & %s " % (dataString(t0,count0),dataString(t1,count1))),

        
    print("\\\\")
print "\\hline"
print "\\end{tabular}"
