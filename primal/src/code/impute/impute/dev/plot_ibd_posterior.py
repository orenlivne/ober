#!/usr/bin/env python
'''
============================================================
Plot our HMM IBD posterior probabilities for model cases. 

Created on December 19, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''

import matplotlib.pyplot as P, os, numpy as np
from numpy.core.function_base import linspace
from impute.dev import ibd_single

N = 50
p = linspace(0, 1, N) + 1e-15
q = 1 - p
colors = ['b', 'g', 'r', 'm', 'k']

def plot_shared_allele_inflation(f):
    '''Plot the phased-phased IBD posterior probability change when we use shared allele frequencies.'''
    P.clf()
    for f in [0.1, 0.5, 0.9]:
        OR = f / (1 - f)
        ibd1 = 1 / (1 + p / OR)
        ibd2 = 1 / (1 + (p ** 2 + q ** 2) / OR)
        P.plot(p, ibd1 / ibd2, label='$f=%.1f,\, OR=%.1f$' % (f, OR))
    P.grid(True)
    P.xlim([0, 1])
    P.ylim([0, 2])
    P.xlabel('Shared Allele Frequency ($p$)')
    P.ylabel(r'$\frac{P\left(IBD|IBS,H_1=h\right)}{P\left(IBD|IBS\right)}$', fontsize=16)
    P.legend()
    P.title('IBD Posterior Probability Inflation Due to Shared Allele')
    
def plot_ibd_posterior(f, Delta_function):
    '''Compare phased-phased and unphased-unphased IBD posteriors. Delta is a kinship model functor.'''
    ibd_hh = ibd_single.ibd_posterior_hh(f, p)
    T = [1, 3, 5, 7]
    j = [(0, 0), (0, 1), (1, 0), (1, 1)]
    ibd_gg = [ibd_single.ibd_posterior_gg(np.array(Delta_function(f)), p, t) for t in T]

    P.clf()
    P.hold(True)
    P.plot(p, ibd_hh, colors[0] + '-', label='HH')
    for x in xrange(len(T)):
        jx = j[x]
        P.plot(p, ibd_gg[x], colors[x + 1] + '--', label='GG, $(%d,%d)$' % (jx[0], jx[1]))
    #P.plot(p, ibd_pu, c + '-,')
    P.grid(True)
    P.xlim([0, 1])
    P.ylim([0, 1])
    P.xlabel('Shared Allele Frequency ($p$)')
    P.ylabel(r'$P\left(IBD|IBS\right)$', fontsize=16)
    P.legend(loc='lower left')
    P.title('IBD Posterior Probability, $f=%.3f$' % (f,))
    return ibd_hh, ibd_gg

#---------------------------------------------
# Main Program
#---------------------------------------------           
if __name__ == '__main__':
    '''Main program - accepts CLI arguments.'''

    out_dir = os.environ['OBER'] + '/doc/ibd'
    P.figure(1)
    plot_shared_allele_inflation([0.1, 0.5, 0.9])
    P.savefig('%s/shared_allele_change.eps' % (out_dir,))

    models = {
               'outbred': lambda f: [0, 0, 0, 0, 0, 0, 0, 4 * f, 1 - 4 * f],
               #'uniform': lambda f: [0.2 * (1 - f), 0.25 * f, 0.2 * (1 - f) , 0.25 * f, 0.2 * (1 - f), 0.25 * f, 0.2 * (1 - f), 0.2 * (1 - f), 0.25 * f]
               'inbred': lambda f: [0.04, 0.04, 0.02, 0.06, 0.04, 0.04, 0.04, 4 * f - 3.5 * 0.04, 1 - 4 * f - 3.5 * 0.04]
               }
    
    fig_num = 2
    for (kinship_model, Delta) in models.iteritems():
        for f in [1. / 32, 1. / 8]:
            P.figure(fig_num)
            plot_ibd_posterior(f, Delta)
            P.savefig('%s/ibd_posterior_%s_f=%.3f.eps' % (out_dir, kinship_model, f))
            #P.title('IBD Posterior Probability, %s, $f=%.3f$' % (kinship_model.title(), f))
            fig_num += 1
        
    P.show()
