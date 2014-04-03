#!/usr/bin/env python
'''
============================================================
Plot IBD HMM posterior -- haplotypes.

Mark:
'Well, the raw IBDLD output gives those posterior probabilities
(posterior because it conditions on all the genotype data). If you
only have the output for the identified segments, that just means the
probability for all SNPs in that segment was greater than some
threshold. The threshold is chosen at the time the program was run and
should be recorded somewhere (I think it was probably set to something
like 0.9 or 0.95). So, a rough model would put the endpoints at this
threshold and somewhat above it in the middle. If the threshold was
0.5 the middle could be anything from 0.51 to 1.0. If it's a long
segment, probably close to 1.0 while a short segment might be more
like 0.51.'

Lide Han says he uses 0.9.
 
Created on Jul 12, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, matplotlib.pyplot as P, os
from impute.ibd.distant.ibd_hmm import ProbIbdHmmCalculator
from impute.ibd.distant.ibd_hmm_hap import kinship, ProbIbdHmmHapCalculator

# Different kinship models
models = {
          'sibs': (0.3, np.array([0, 0, 0, 0, 0, 0, 0.25, 0.5, 0.25]))
#          'uniform':  (0.3, np.ones(9) / 9.0)
          }

def ibd_posterior_gg(lam, Delta, m):
    '''Simple case: sibs with (1,1) genotypes along the segment.'''
    x = np.arange(m)
    m = len(x)
    return ProbIbdHmmCalculator(lam=lam, Delta=Delta, x=x, p=0.2 * np.ones(m),
                                g1=np.ones((m, 2), dtype=np.uint),
                                g2=np.ones((m, 2), dtype=np.uint),
                                e=0.01).prob()

def ibd_posterior_hh(lam, Delta, m):
    '''Simple case: sibs with (1,1) genotypes along the segment.'''
    x = np.arange(m)
    m = len(x)
    return ProbIbdHmmHapCalculator(lam=lam, f=kinship(Delta), x=x, p=0.2 * np.ones(m),
                                   h1=np.zeros((m,), dtype=np.uint),
                                   h2=np.zeros((m,), dtype=np.uint),
                                   e=0.01).prob()
    
def plot_posterior(lam, Delta, n, i, m):
    '''Create a confidence vs. x plot for a particular segment size M.'''
    gg = ibd_posterior_gg(lam, Delta, m)
    hh = ibd_posterior_hh(lam, Delta, m)
    P.subplot(n, 1, i)
    x = np.arange(m)
    P.plot(x, gg, 'b', x, hh, 'r')
    P.ylabel('P(IBD)')
    P.ylim(0.7, 1.0)
    P.xlim(x[0], x[-1])
    P.xlabel('m = %.1f' % (m,))
    #print '#SNPs m = %.2f' % (m,)
    return (gg, hh)

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    out_dir = os.environ['OBER'] + '/doc/ibd/hmm-hap'
    
    np.set_printoptions(precision=15, suppress=True)    
    m = 1024
    lam, Delta = models['sibs']
    gg = ibd_posterior_gg(lam, Delta, m)
    hh = ibd_posterior_hh(lam, Delta, m)
    a = [gg, hh]
    print m, np.max(gg), np.max(hh)
    
    # Create a confidence vs. x plot for various segment sizes M
    np.set_printoptions(precision=2, suppress=True)
    n = 5
    fig_num = 0
    for model, (lam, Delta) in models.iteritems():
        fig_num += 1
        P.figure(fig_num)
        P.clf()
        [plot_posterior(lam, Delta, n, i, 4 ** i) for i in xrange(1, n + 1)] 
        
        P.subplot(n, 1, 1)
        P.title('IBD Posterior Probability\n$\lambda = %.2f, \Delta = %s, f = %.2f$' % (lam, repr(Delta)[6:-1], kinship(Delta)))
        P.subplot(n, 1, n)
        
        P.xlabel('Genetic Distance')
        # Crappy aspect ratio but doesn't matter for now
        P.show() 
        P.savefig('%s/ibd_%s.png' % (out_dir, model,))
