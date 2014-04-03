#!/usr/bin/env python
'''
============================================================
Validate our detailed identity coefficients by summing them
up to condensed coefficients and comparing with Lide's
IBDLD results.

Created on March 17, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, os, matplotlib.pyplot as P  # , db_gene
from scipy import stats
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#---------------------------------------------
# Constants
#---------------------------------------------
'''Number of samples.'''
n = im.hutt_pedigree().num_genotyped

'''Converts detailed to condensed by multiplying by it on the right.'''
A = np.matrix(np.zeros((15, 9)))
S = im.poo.idcoef.S
A[np.arange(len(S)), S] = 1

#---------------------------------------------
# Methods
#---------------------------------------------
def row(i, j):
    '''Delta array row corresponding to the sample pair (i,j).'''
    return i * n + j

def ind(r):
    '''Sample pair (i,j) of row r.'''
    return divmod(r, n)

def error_data(ibdld, detailed):
    '''Return the error in delta7 and the kinship f between IBDLD and PRIMAL''s detailed coefficients.'''
    ibdld_delta7, ibdld_f = ibdld[:, 0], ibdld[:, 1]
    delta = detailed * A
    f = im.ibd_hmm_hap.kinship_vector(delta)
    return ibdld_delta7 - delta[6], ibdld_f - f

def primal_data(detailed):
    '''Return delta7, kinship f from PRIMAL''s detailed coefficients.'''
    delta = detailed * A
    return np.array(np.concatenate((delta[:, 6], im.ibd_hmm_hap.kinship_vector(delta)), axis=1)) 

def plot_comparison(x, y, fig_num=1, max_y=1, xlabel='x', ylabel='y', title='Comparison'):
    '''Scatter-plot the two corresponding statistics x and y produced by IBDLD and PRIMAL, respectively.'''   
    a, b, r, _, stderr = stats.linregress(x, y)
    xr = np.polyval([a, b], x)
    P.figure(fig_num)
    P.clf()
    P.hold(True)
    P.grid(True)
    P.scatter(x, y, lw=0, s=2, color='b')
    P.plot(x, xr, 'r-')
    X = np.arange(0, max_y, 0.01)
    P.plot(X, X, 'g-')
    P.xlim([-0.01, max_y])
    P.ylim([-0.01, max_y])
    P.xlabel(xlabel)
    P.ylabel(ylabel);
    P.title('%s\nslope = %.2f, intercept = %.2f, r = %.2f, std error = %.2e' % (title, a, b, r, stderr))

####################################################################################
if __name__ == '__main__':
    # Load Delta7, kinship (phi) from Lide Han's IBDLD v3 output (see the compare-segments-v3 script)
    ibdld = np.loadtxt(os.environ['OBER_OUT'] + '/ibd/ibdld-v3/ibdld.out', dtype=np.float32)
    
    # Load our condensed coefficients
    # ibd = im.index.segment_index.SegmentIndex(os.environ['OBER_OUT'] + '/index_segments')
    # t_counter, total_bp = im.poo.idcoef.idcoef(ibd, 22, algorithm='full', out='count', do_region=0)
    # t0 = im.poo.idcoef.reshape_idcoef_array(t_counter)
    # detailed = (t0.astype(float) / total_bp) * A
    detailed = np.loadtxt(os.environ['OBER_OUT'] + '/po/idcoefs-all.txt', dtype=np.float32)
    primal = primal_data(detailed)
    
    # Compare the two methods
    E = primal - ibdld
    # f, E = error_data(ibdld, detailed)
    plot_comparison(ibdld[:, 1], primal[:, 1], fig_num=1, max_y=0.6, xlabel='IBDLD Kinship $\phi$', ylabel='PRIMAL Kinship \phi', title='PRIMAL-IBDLD Kinship Comparison')
    P.savefig(os.environ['OBER'] + '/doc/poo/primal-ibdld-comparison-kinship.png')
    plot_comparison(ibdld[:, 0], primal[:, 0], fig_num=2, max_y=1.05, xlabel='IBDLD $\Delta_7$', ylabel='PRIMAL $Delta_7$', title='PRIMAL-IBDLD $\Delta_7$ (IBD=2 State) Comparison')
    P.savefig(os.environ['OBER'] + '/doc/poo/primal-ibdld-comparison-delta7.png')
    # plot_error(f, E)
