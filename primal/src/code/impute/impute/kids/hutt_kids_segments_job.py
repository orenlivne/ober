#!/usr/bin/env python
'''
============================================================
A script for calculating IBD segments of one Hutt kid sib
vs. all other affy individuals. To be parallelized on a
cluster.

Created on August 5, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, numpy as np, sys
from impute.kids.hutt_kids import ibd_segments, get_sib_ids, read_chip_problem

####################################################################################
if __name__ == '__main__':
    chip, sib = sys.argv[1:]

    path = os.environ['OBER_OUT'] + '/kids'
    chrom = 22
    debug = 1
    
    sibs = get_sib_ids(path)
    problem = read_chip_problem(path, chrom, chip)    
    affy_samples = np.array(list(set(xrange(problem.pedigree.num_genotyped)) - set(sibs)))
    ibd_segments(problem, sib, affy_samples, path + '/%s/chr%d/segments_%d.out' % (chip, chrom, sib))
