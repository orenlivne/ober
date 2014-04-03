#!/usr/bin/env python
'''
============================================================
IMPUTE2 QC: plot het concordance vs. IMPUTE2 probability
threshold. 

Created on February 14, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, matplotlib.pylab as P, numpy as np

def plot_het_concordance(a):    
    x = P.linspace(0, 1, a.shape[0])
    y = a[:, 1] / (a[:, 0] + 1e-15)
    P.figure(1)
    P.clf()
    ax = P.subplot(111, title='IMPUTE2-PRIMAL Discordance Rate', xlabel='IMPUTE2 Confidence', ylabel='Discordance Rate')
    ax.plot(x, y, 'b.-')
    P.xlim([min(x[-20:]), 1.01])
    P.grid(True)
    P.ylim([0, max(y[-20:])])
    P.hold(True)
    
    f = np.polyfit(1 - x[-10:], y[-10:], 1)
    yFit = np.polyval(f, 1 - x)    
    # ax.plot(x, 0.3 * (1 - x) + 0.005, 'g--')  # Regression line
    ax.plot(x, yFit, 'g--')  # Regression line
    
    for item in [ax.title, ax.xaxis.label, ax.yaxis.label]: item.set_fontsize(18)
    for item in (ax.get_xticklabels() + ax.get_yticklabels()): item.set_fontsize(14)
    return x, y, f
    
####################################################################################
if __name__ == '__main__':
    x, y, f = plot_het_concordance(np.loadtxt(open(sys.argv[1], 'rb'), dtype=int))
    print 'Regression line: het.concord ~ %.4f (1-prob) + %.4f' % (f[0], f[1])
    if len(sys.argv) >= 3:
        P.savefig(os.environ['OBER'] + '/doc/qc-impute2/' + sys.argv[2])
