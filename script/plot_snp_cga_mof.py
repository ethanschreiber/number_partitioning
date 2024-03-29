from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import ProcessExperimentsLibrary  # Import Library
lib = ProcessExperimentsLibrary   # Give short name
from optparse import OptionParser

import matplotlib.pyplot as plt
import numpy as np

def plotExperiments(filename, MARKERS,ALG_NAMES, numBits, minN, maxN,  ALGORITHMS, COLOR): 
    minK=3
    maxK=12
    outFile = open(filename,"w")
    
    NArray=[]
    for N in range (minN,maxN+1) :
        NArray.append(N)
    
    fig, axarr = plt.subplots(5,2, sharex=True)
    fig.set_canvas(plt.gcf().canvas)
    fig.set_size_inches(7.5,9)
    coordsX = 0
    coordsY = 0

    legend = []
    for i in range(len(ALGORITHMS)) :
        legend.append(ALG_NAMES[ALGORITHMS[i]])
    for k in range(minK,maxK+1) :
        ax = axarr[coordsX,coordsY]

        allTimes = []
        for i in range(len(ALGORITHMS)) :
            allTimes.append([])
            
        for N in range(minN, maxN+1) :
            for i in range(len(ALGORITHMS)) :
                algorithm = ALGORITHMS[i]
                inputFilename = lib.filename(k,N,algorithm,100,numBits)
                print(inputFilename)
                time = lib.processExperiment(inputFilename,100,1)
                allTimes[i].append(time)

                    
        ax.set_title('k = ' + str(k))
                    
        # See http://matplotlib.org/1.3.1/api/axes_api.html#matplotlib.axes.Axes.plot
        # For styles

        for i in range(len(ALGORITHMS)) :
            algorithm = ALGORITHMS[i]
            marker = MARKERS[algorithm]            
            color = COLOR[algorithm]
            #plt.plot(NArray,allTimes[i],marker,markerfacecolor='w',markersize=5)
            #axarr[k-minK].plot(NArray,allTimes[i])
            ax.plot(NArray,allTimes[i],marker,markerfacecolor=color,markersize=5)
        
        if k == minK :
            ax.legend(legend, loc='upper left', fontsize=10)
        if (coordsX == 4) :        
            ax.set_xlabel("n")
        
        if (coordsY == 0) :
            ax.set_ylabel("Time (s)")


        # Make labels smaller
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.tick_params(axis='both', which='minor', labelsize=8)

        ax.set_yscale('log')
        ax.grid(True) 
        #maxValue = max(map(max,allTimes))
        #minValue = 0 - (maxValue / 50.0)
        
        #plt.ylim(minValue,maxValue)

        if coordsY == 1 :
            coordsX = coordsX + 1
            coordsY = 0
        else :
            coordsY = 1

            #plt.suptitle('Average Run Time for 48-Bit Instances')
    print("Writing",filename)    
        #pdf.savefig(bbox_inches='tight')

    fig.savefig(filename, format='pdf',bbox_inches='tight')
    
    
# =====================================
# Run Program for cga, snpess and snpie
# =====================================

filename = "/home/ethan/workspace/journal/images/snp_cga_moffitt_3_10_plot.pdf"

MARKERS={"cga_mw" : "ks--",\
         "snp" : "ko--",\
         "moffitt"  : "k^--"}

ALG_NAMES = {"cga_mw" : "CGA",\
             "snp" : "SNPESS",\
             "moffitt"  : "SNPIE"}

COLOR = {"cga_mw"     : '1.0',\
         "snp"     : '0.5',\
         "moffitt" : '0.0'}

plotExperiments(filename,MARKERS,ALG_NAMES,48,20,45,['cga_mw','moffitt','snp'],COLOR)
         
# ===================================================
# Run Program for cga, snpess, snpie, BSIBC and BSBCP
# ===================================================

filename = "/home/ethan/workspace/journal/images/snp_moffitt_bsbcp_bsibc_3_10_plot.pdf"

MARKERS = {"snp"      : "ks--",\
           "moffitt"  : "ks--",\
           "bsbc_IE50": "ko--",\
           "bsbcp"    : "ko--"}           

ALG_NAMES = {"snp"      : "SNPESS",\
             "moffitt"  : "SNPIE",\
             "bsbc_IE50": "BSIBC",\
             "bsbcp"    : "BSBCP"}

COLOR = {"snp"      : '0.5',\
         "moffitt"  : '0.0',\
         "bsbc_IE50": '1.0',\
         "bsbcp"    : '0.0'}    

plotExperiments(filename,MARKERS,ALG_NAMES,48,20,45,['moffitt','snp','bsbc_IE50','bsbcp'],COLOR)
