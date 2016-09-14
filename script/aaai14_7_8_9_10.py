from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os.path
import matplotlib.pyplot as plt
import numpy as np



def plot(filename,k,X,alg1,alg2,name1,name2) :
    plt.rcParams.update({'font.size': 25})
    pdf = PdfPages(filename) 
    plt.figure(figsize=(10,7))
    plt.subplot(1,1,1)
    plt.yscale('log')
    
    print (len(X))
    print (len(alg1))
    print (len(alg2))
    plt.plot(X,alg1)
    plt.plot(X,alg2)
    

    plt.axis([X[0], X[-1], 0, alg2[-1]])
    
    plt.grid(True)
    plt.xlabel('n')
    plt.ylabel('time (s)')
    plt.legend([name1,name2], loc='upper left')

    plt.title("k=" + str(k),fontsize=40)
    print("Writing",filename)    
    pdf.savefig()
    plt.close()
    pdf.close()
    plt.show()


X = range(40,61)

CIW07=[.219,.207,.292,.407,.522,.753,.899,1.41,1.61,2.67,2.93,5.43,6.36,10.6,12.2,24.2,27.1,46.6,53.0,106,119]
SNP07=[1.96,3.02,5.03,7.47,11.1,19.5,26.7,41.4,62.9,92.9,135,233,364,529,864,1328,1945,3158,5015,7520,10715]
CIW08=[.189,.226,.330,.433,.510,.677,.820,1.25,1.48,2.20,2.29,4.25,5.08,8.45,9.33,18.4,20.6,34.7,41.1,71.2,80.9]
MOF08=[1.24,2.24,3.77,5.88,8.45,14.3,23.7,38.2,60.0,100,166,281,418,707,1189,2082,3078,5589,9182,13441,23085]
CIW09=[.275,.247,.344,.447,.506,.713,.935,1.39,1.55,2.04,2.23,3.42,3.97,6.77,8.12,14.5,16.7,23.8,29.1,48.9,56.6]
MOF09=[1.25,2.08,3.26,5.35,8.65,13.7,24.8,41.4,59.4,98.0,154,274,382,724,1116,1701,2878,4983,7636,11874,21414]
CIW10=[.404,.441,.630,.644,.840,.968,.913,1.35,1.42,1.93,2.32,3.82,5.28,7.56,8.02,14.5,14.4,18.0,23.7,35.6,42.0]
MOF10=[1.39,2.03,3.27,5.27,9.25,14.0,21.9,42.6,56.5,92.8,141,263,362,735,1120,1920,2700,4950,8508,12249,18036]

DIR = '/home/ethan/workspace/aaai14/presentation/'

plot(DIR + "ciw07.pdf",7 ,X,CIW07,SNP07,'CIW (AAAI-14)','SNP (ICAPS-13)')
plot(DIR + "ciw08.pdf",8 ,X,CIW08,MOF08,'CIW (AAAI-14)','MOF (IJCAI-13)')
plot(DIR + "ciw09.pdf",9 ,X,CIW09,MOF09,'CIW (AAAI-14)','MOF (IJCAI-13)')
plot(DIR + "ciw10.pdf",10,X,CIW10,MOF10,'CIW (AAAI-14)','MOF (IJCAI-13)')

