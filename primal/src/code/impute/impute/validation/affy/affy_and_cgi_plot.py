#!/usr/bin/env python
'''
============================================================
Imputation performance plot for Carole's Hutterites grant
renewal.

Created on Oct 19, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, os, matplotlib.pylab as P

a = np.load(os.environ['OBER_OUT'] + '/validation/call_rate.npz')
x = a['maf_bins']
maf = np.array([0.5 * (x[i] + x[i + 1]) for i in xrange(len(x) - 1)])

P.figure(1)
P.clf()
P.grid(True)
P.hold(True)
P.plot(a['stats_maf'], a['stats_concordance'], 'r-', label='Concordance (all genotypes)')
P.plot(a['stats_maf'], a['stats_het'], 'g-', label='Concordance (het genotypes)')
P.plot(maf, a['mean_call_rate'][1:], 'b-', label='Call Rate')
P.legend(loc='lower right', prop={'size': 12})
P.xlabel('Minor Allele Frequency')
P.title('Imputation Performance')

P.savefig(os.environ['OBER'] + '/doc/imputation/cgi/performance.png')
