#!/usr/bin/env python
'''
============================================================
Plot phasing % after the different phasing stages. 

Created on August 16, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, matplotlib.pyplot as P, numpy as np, os
from numpy.core.function_base import linspace

P.figure(1)
P.clf()
P.hold(True)

p = im.hutt('hutt.phased.npz')
d5 = im.plots.plot_fill_fraction(p, color='b', label='Stage 5')

p = im.hutt('hutt.stage6.npz')
d = im.plots.plot_fill_fraction(p, color='r', label='Stage 6')

zoom = 0.96
ticks = 10
min_y = 0.95 * min(d[:, 1][0], d5[:, 1][0])
max_x = np.where(d[:, 1] > zoom)[0][0]
P.xlim([0, max_x + 1])
P.ylim([min_y, 1.0])
yticks = linspace(min_y, 1.0, ticks)
P.yticks(yticks, ['%.3f' % (t,) for t in yticks])

P.title('Hutterites Phasing Coverage, Chromosome 22')
P.legend(loc='lower right', prop={'size': 12})
P.show()

P.savefig(os.environ['OBER'] + '/doc/phasing/fill.png')
