from __future__ import print_function
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os.path
import matplotlib.pyplot as plt
import numpy as np


def plot(filename,X,RNP,IRNP,MOF,CIW,plots) :
    pdf = PdfPages(filename) 
    plt.figure(figsize=(8,6))
    plt.subplot(1,1,1)
    plt.yscale('log')
    
    plt.plot(X, RNP, 'r')
    plt.text(25.5,5000,'RNP\nKorf, IJCAI-09',color='r')

    if plots >= 2 :
        plt.plot(X, IRNP, 'g')
        plt.text(36,5000,'IRNP\nKorf, IJCAI-11',color='g')
    else :
        plt.plot(X, IRNP, 'w',markeredgecolor='w')  # for axes to show correctly
        
    if plots  >= 3 :
        plt.plot(X, MOF, 'y')
        plt.text(49,5000,'MOF\nMoffit, IJCAI-13',color='y')

    if plots >= 4 :
        plt.plot(X, CIW, 'b')
        plt.text(49,.5,'CIW\nSchreiber, AAAI-14',color='b')

    plt.axis([20, X[-1], 0, 18500])
    
    plt.grid(True)
    plt.xlabel('n',fontsize=20)
    plt.ylabel('time (s)',fontsize=20)
    #plt.legend(['RNP (IJCAI 09)', 'IRNP (IJCAI 11)', 'MOF (IJCAI 13)', 'CIW (AAAI 14)'], loc='upper left')
    plt.title('Time to Optimally Partition n \nIntegers into k=10 Subsets',fontsize=20)
    print("Writing",filename)    
    pdf.savefig()
    plt.close()
    pdf.close()
    plt.show()

# evenly sampled time at 200ms intervals

RNP  = [57, 87, 385, 590, 1388, 1971]
IRNP = [0 ,0 ,.001,.011,.041,.061,1,2,8,28,70,693,917,2407,4079,6837]
MOF  = [.000,.000,.000,.000,.000,.000,.001,.001,.003,.006,.012,.062,.022,.051,.099,.157,.232,.531,.747,1,1.39,2.03,3.27,5.27,9.25,14.0,21.9,42.6,56.5,92.8,141,263,362,735,1120,1920,2700,4950,8508,12249,18036]
CIW  = [.004,.017,.027,.055,.158,.246,.408,3.00,2.14,.536,1.77,1.1,.168,.119,.150,.104,.112,.115,.182,.174,.404,.441,.630,.644,.840,.968,.913,1.35,1.42,1.93,2.32,3.82,5.28,7.56,8.02,14.5,14.4,18.0,23.7,35.6,42.0]
        
while (len(RNP) < len(MOF)) :
    RNP.append(100000000)

while (len(IRNP) < len(MOF)) :
    IRNP.append(100000000)


X = range(20,20+len(RNP))

DIR = '/home/ethan/workspace/aaai14/presentation/'
plot(DIR + "RNP.pdf"             ,X,RNP,IRNP,MOF,CIW,1)
plot(DIR + "RNP_IRNP.pdf"        ,X,RNP,IRNP,MOF,CIW,2)
plot(DIR + "RNP_IRNP_MOF.pdf"    ,X,RNP,IRNP,MOF,CIW,3)
plot(DIR + "RNP_IRNP_MOF_CIW.pdf",X,RNP,IRNP,MOF,CIW,4)
