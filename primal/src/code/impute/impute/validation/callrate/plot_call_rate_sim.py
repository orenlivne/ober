#!/usr/bin/env python
'''
============================================================
Plotting tool accompaniment to call_rate_sim.py

Created on March 12, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib.pyplot as P, numpy as np

####################################################################################
if __name__ == '__main__':
    a = np.loadtxt('/home/oren/ober/out/pedigree/call-rates-10000.txt')
    i = np.argsort(a[:, 1])
    P.figure(1)
    P.clf()
    P.plot(a[i, 1], 'g', label='Genotype')
    P.plot(a[i, 2], 'b', label='Allele') 
    P.title('Ideal Hutterites Call Rates')
    P.xlabel('Sample #') 
    P.ylabel('Call Rate')
    x = P.xlim()
    P.xlim([-10, x[1]])
    P.legend(loc='lower right')
    P.savefig('call-rate-sample.png')
    
    b = np.mean(a, axis=0)
    al = np.linspace(0.9, 1, (a.shape[1] - 1) / 2)
    P.figure(2)
    P.clf()
    P.plot(al, b[1::2], 'go-', label='Genotype')
    P.plot(al, b[2::2], 'bo-', label='Allele')
    P.legend(loc='lower right')
    P.ylim([0.8, 1])
    P.xlabel('Phasing %')
    P.ylabel('Average Sample Call Rate')
    P.title('Average Call Rate vs. Phasing %')
    P.savefig('call-rate-phasing.png')
