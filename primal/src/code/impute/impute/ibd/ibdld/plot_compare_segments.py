#!/usr/bin/env python
'''
============================================================
Plot segment comparison statistics produced by
compare_segments.py for chromosome 22.

Created on February 14, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, matplotlib.pyplot as P, os, sys
from scipy import stats
from numpy.lib.polynomial import polyval
from db_gene.snp.file_dao import ChromDao
# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)

def plot_ibd_total(c):
    '''Plot Oren's total genomic IBD % vs. kinship, displaying the trend of IBD.'''     
    i = np.where(c[:, 4] > 0)[0]
    t, ibd = c[i, 3], (1.0 * c[i, 4]) / sum(ChromDao.TOTAL_BP_TYPED)
    # Fix HBD probabilities to account for the probability that one can choose the same haplotype twice.
    # (We only account for the probability of the two different haps being equal in c[i,4].)
    hbd = np.where(c[:, 0] == c[:, 1])[0]
    ibd[hbd] = 0.5 * (1 + ibd[hbd])
    
    j = i[np.where(c[:, 0] != c[:, 1])[0]]
    a, b, r, _, stderr = stats.linregress(t[j], ibd[j])
    xr = polyval([a, b], t[j])
    
    P.figure(1)
    P.clf()
    P.hold(True)
    P.scatter(t, ibd, s=2, lw=0, color='b')
    P.plot(t[j], xr, 'r-')
    T = np.arange(0, 1, 0.01)
    P.plot(T, T, 'g-')
    P.xlim([-0.01, 1.01])
    P.ylim([-0.01, 1.01])
    P.title(' %% IBD of Whole Genome vs. Expected %%\nslope = %.2f, intercept = %.2f, r = %.2f, std error = %.2e' % (a, b, r, stderr))
    P.xlabel('$\gamma := \Delta_1+\Delta_3+\Delta_5+\Delta_7+\Delta_8$: Empirical ID Coefficients')
    P.ylabel('$i = $ % IBD in Oren''s Code')
    P.show()
    return t, ibd

def plot_comparison(c, save_location=None):     
    i = np.where(c[:, 5] > 0)
    f, r = c[i, 2], (1.0 * c[i, 3]) / c[i, 5]
    P.figure(1)
    P.clf()
    P.scatter(f, r, s=2, lw=0, color='b', label='% IBDLD covered')
    P.xlim([-0.01, 0.6])
    P.ylim([-0.01, 1.01])
    P.xlabel('Kinship')
    P.ylabel('% Segments Covered')
    P.title('% of IBDLD Segments Covered by Oren, Chromosome 22')
    P.show()
    if save_location:
        P.savefig(save_location + '/ibdld-covered.png')  

    i = np.where(c[:, 6] > 0)
    f, r = c[i, 2], (1.0 * c[i, 3]) / c[i, 6]
    P.figure(2)
    P.clf()
    P.scatter(f, r, s=2, lw=0, color='b', label='% Oren covered')
    P.xlim([-0.01, 0.6])
    P.ylim([-0.01, 1.01])
    P.xlabel('Kinship')
    P.ylabel('% Segments Covered')
    P.title('% of Oren Segments Covered by IBDLD, Chromosome 22')
    P.show()
    if save_location:
        P.savefig(save_location + '/oren-covered.png')  

#####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Chromsome 22 data
    # File row format: id1 id2 kinship |A & B| |A | B| |A| |B| where A=IBDLD segments, B=our segments 
#    c = np.loadtxt(os.environ['OBER_OUT'] + '/ibd/comparison.chr22')
#   plot_ibd_trend(c)
#    c = np.loadtxt(os.environ['OBER_OUT'] + '/ibd/total_ibd.out')
#    plot_ibd_total(c)
#    P.savefig(os.environ['OBER'] + '/doc/ibd/ibd-total.png')  

    c = np.loadtxt(sys.argv[1]) #np.loadtxt(os.environ['OBER_OUT'] + '/ibd/total_ibd_empirical.out')
    t, ibd = plot_ibd_total(c)
    P.savefig(os.environ['OBER'] + '/doc/ibd/ibd-total-empirical.png')  
    