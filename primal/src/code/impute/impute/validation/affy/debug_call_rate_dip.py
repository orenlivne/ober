#!/usr/bin/env python
'''
============================================================
Investigating call rate dips far from chromosome edges.
Are they near centromeres?

Created on December 9, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os, numpy as np
from scipy import ndimage
from impute.plot.plots import genome_bp_offset

####################################################################################
if __name__ == '__main__':

    proband, parent = 579, 1403  # 248, 547#80, 547
    
    data = os.environ['OBER_OUT'] + '/validation/affy'
    stats = im.imputation.impute_stats.Stats.load(data + '/stats.npz')
    
    chrom = stats.snp.info['chrom']
    x = stats.cm_edge_dist
    y = stats.cm_cumulative
    c = stats.snp.call_rate_imputed_full
    k = np.argsort(x)
    cs = ndimage.median_filter(c, 20)
    dip = np.where((cs[k] < 0.58) & (x[k] > 70))[0]
    kd = k[dip]
    i = np.argsort(chrom[k[dip]])
    genome_offset = genome_bp_offset()[:22]
    for z in zip(chrom[kd][i], stats.snp.info['base_pair'][kd][i], stats.snp.info['base_pair'][kd][i] + genome_offset[chrom[kd][i] - 1], x[kd][i], y[kd][i]):
        print '%-2d  %-10d %-12d  %-9.3f %-9.3f' % (z)
