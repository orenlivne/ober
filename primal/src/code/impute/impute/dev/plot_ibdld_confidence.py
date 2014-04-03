#!/usr/bin/env python
'''
============================================================
IBDLD confidence model = approximate posterior HMM
probabilities. Based on discussion with Mark Abney.

Mark:
'Well, the raw IBDLD output gives those posterior probabilities
(posterior because it conditions on all the genotype data). If you
only have the output for the identified segments, that just means the
probability for all SNPs in that segment was greater than some
threshold. The threshold is chosen at the time the program was run and
should be recorded somewhere (I think it was probably set to something
like 0.9 or 0.95). So, a rough model would put the end points at this
threshold and somewhat above it in the middle. If the threshold was
0.5 the middle could be anything from 0.51 to 1.0. If it's a long
segment, probably close to 1.0 while a short segment might be more
like 0.51.'

Lide Han says he uses 0.9.
 
Created on Jul 12, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''

import numpy as np
import matplotlib.pyplot as plt
from impute.ibd.ibdld import ibd_ld
    
def plot_ibd_ld(n, i, M):
    '''Create a confidence vs. x plot for a particular segment size M.'''
    x = np.linspace(-0.5 * M, 0.5 * M, 100)
    plt.subplot(n, 1, i)
    plt.plot(x, ibd_ld.ibd_ld_confidence(x, M))
    plt.ylabel('IBD Probability')
    plt.ylim(ibd_ld.IBDLD_THRESHOLD, 1.0)
    plt.xlabel('M = %.1f' % (M,))
    print 'M = %.2f' % (M,)

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Create a confidence vs. x plot for various segment sizes M
    plt.clf()
    plt.figure(1)
    n = 4
    [plot_ibd_ld(n, i, 10 ** (i - 2)) for i in xrange(1, n + 1)] 
    
    plt.subplot(n, 1, 1)
    plt.title('IBD Probability Model for Various Segment Lengths')
    plt.subplot(n, 1, n)
    
    plt.xlabel('SNP position')
    # Crappy aspect ratio but doesn't matter for now
    plt.show() 
    plt.savefig('ibd_confidence.png')
    