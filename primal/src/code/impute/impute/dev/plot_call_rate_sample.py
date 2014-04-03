#!/usr/bin/env python
'''
============================================================
Plotting tool accompaniment to call_rate_sim.py

Created on March 12, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib.pyplot as P, numpy as np, os

####################################################################################
if __name__ == '__main__':
    a = np.loadtxt(os.environ['OBER'] + '/doc/imputation/cgi/v2/sample.txt')
    P.figure(1)
    P.clf()
    P.plot(np.sort(a[:, 2]))
    P.xlabel('Sample #')
    P.ylabel('Call Rate')
    P.title('Hutterites CGI Imputation Call Rate')
    P.savefig(os.environ['OBER'] + '/doc/imputation/cgi/v2/call-rate-sample.png')
