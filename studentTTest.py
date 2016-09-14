import csv
import sys
import string
    
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
print "\n\n"

CSV = [] # Matrix to store values

FILENAMES = sys.argv[1:] # The input files start as the second argument, program name is first

OUTPUT_FILENAME="./studentTTest.csv"
outFile = open(OUTPUT_FILENAME,'w')



# Read from files and write to CSV matrix
for fileIdx in range(0,len(FILENAMES)) :
    CSV.append([])  # Add row
    print "File Idx:", fileIdx
    CSV[fileIdx].append(FILENAMES[fileIdx])
    print "Processing file: ", FILENAMES[fileIdx]
    inFile  = open(FILENAMES[fileIdx], "rb")
    reader = csv.reader(inFile)
        
    for row in reader:    
        print "row: ", row
        print row[1]
        print "File Idx:", fileIdx        
        print len(CSV)
        CSV[fileIdx].append(row[1])        
    
    inFile.close()

rowNum = 0
foundOne = True


while foundOne :
    foundOne = False
    for i in range (0, len(CSV)) :
        if i > 0 :
            outFile.write(", ")
        if rowNum < len(CSV[i]) :
            outFile.write(CSV[i][rowNum].rjust(10))
            foundOne = True
        else:
            outFile.write("".rjust(10))

    outFile.write("\n")
    rowNum += 1

print "Wrote file: ",OUTPUT_FILENAME
