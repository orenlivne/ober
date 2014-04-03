#!/usr/bin/env python
'''
============================================================
Extract imputed Hutterite genotypes from a tabixed imputed
genotype file set.

Prerequisites: tabix

Created on April 17, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
# Compare imputation call rate for Minal's SNPs before and after LD-based imputation was integrated
import numpy as np, matplotlib.pyplot as P, os

def digitize_into_bins(data, bins):
    '''Return an array with averages of ''data'' over each bin in ''bins''.''' 
    maf_bin = np.digitize(data, bins)
    avg = np.array([(1.*np.mean(data[maf_bin == i])) for i in xrange(len(maf_bins))])
    avg[~np.isfinite(avg)] = 0
    return avg

####################################################################################
if __name__ == '__main__':

    # Number of samples
    N = 1415
    # Old imputation data 
    c = np.loadtxt(os.environ['OBER_OUT'] + '/requests/minal_interacting/minal2.frq', usecols=(4, 5), skiprows=1, dtype=float)
    call_rate1 = c[:, 1] / (2 * N)
    # New imputation data 
    d = np.loadtxt(os.environ['OBER_OUT'] + '/requests/minal_response/minal.frq', usecols=(4, 5), skiprows=1, dtype=float)
    call_rate2 = d[:, 1] / (2 * N)
    # Bin call rates by MAF
    maf_n = 40
    maf_bins = np.linspace(0, 0.5, maf_n + 1)
    
    P.figure(1)
    P.clf()
    P.hold(True)
    P.grid(True)

#    fig, ax = P.subplots()    
#    c1 = digitize_into_bins(call_rate1, maf_bins, label='Pedigree Imputation, mean=%.2f' % (np.mean(call_rate1),)) 
#    P.bar(maf_bins, c1, color='b', width=0.5)
#
#    d1 = digitize_into_bins(call_rate1, maf_bins) 
#    P.bar(maf_bins, d1, color='r', width=0.5, label='Super Imputation, mean=%.2f' % (np.mean(call_rate2),))
    P.scatter(c[:, 0], call_rate1, color='b', label='Pedigree Imputation, mean=%.2f' % (np.mean(call_rate1),)); 
    P.scatter(d[:, 0], call_rate2, color='r', label='Super Imputation, mean=%.2f' % (np.mean(call_rate2),)); 
    P.xlabel('Minor Allele Frequency')
    P.ylabel('Call Rate')
    P.title('Imputation Call Rate: %d SNPs' % (max(c.shape[0], d.shape[0]),))
    P.xlim([0, 0.5])
    P.legend(loc='lower right', prop={'size': 12})
    P.show()
    
    P.savefig(os.environ['OBER_OUT'] + '/requests/minal_response/minal_call_rate.png')
