'''
Created on Jan 31, 2013

@author: oren
'''
import sys, matplotlib.pylab as P, numpy as np

if __name__ == '__main__':
    s = np.loadtxt(sys.argv[1])
    P.figure(1) 
    P.clf()  
    P.hist((s[:, 3] - s[:, 2]) / 1e6, 50)
    P.xlabel('Length [Mbp]')
    P.ylabel('Frequency')
    P.title('IBD Segment Length Distribution in the Hutterites')
