#!/usr/bin/env python
'''
============================================================
CGI imputation call rate histograms by SNP and sample.

The input file produced from count.txt via the command

awk '{ a0 = $4 + $6 + 2*$7 + $8 + $10; a1 = $5 + $8 + $9 + $10 + 2*$11; b = a0+a1+1e-15; if (a0 < a1) { maf = a0/b; } else { maf = a1/b; } printf "%f %f\n", maf, $34; }' count.txt > call-rate.txt

Created on Oct 21, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, impute as im, os, numpy as np, matplotlib.pylab as P

def plot_call_rate(c):
    # Histogram
    P.clf()
    P.figure(1)
    P.hist(c[:,1], normed=True)
    P.xlabel('Call Rate')
    P.ylabel('Portion of Variants')
    P.savefig(os.environ['OBER'] + '/doc/imputation/cgi/call_rate.png')

####################################################################################
#if __name__ == '__main__':
#    # Input parameters
#    file_name = sys.argv[1]  # Name of data file with MAF, call rates
#
#    # Load data
#    c = np.loadtxt(file_name, dtype=np.float16)
#
#    # Breakdown by call rate (proportional to the #samples, 1415)
#    plot_call_rate(c)
#    h = np.histogram(c[:,1])
#    a = np.flipud(np.cumsum(np.flipud(h[0])))/float(c.shape[0])
#    print np.concatenate((h[1][:-1][newaxis].transpose(), a[newaxis].transpose()), axis=1)

    # Breakdown by minor allele frequency
    maf_n = 20
    maf_bins = np.linspace(0, 0.5, maf_n + 1)
    maf_bin = np.digitize(c[:,0], maf_bins)
    d = c.astype(float64)
    mean_call_rate = np.array([(1.*np.mean(d[maf_bin == i,1])) for i in xrange(len(maf_bins))])
    P.bar(maf_bins - h, mean_call_rate, width=h)

    P.figure(2)
    h = (maf_bins[-1] - maf_bins[0]) / maf_n
    P.bar(maf_bins - h, mean_call_rate, width=h)
    P.savefig(os.environ['OBER'] + '/doc/imputation/cgi/call_rate_maf.png')
