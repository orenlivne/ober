#!/usr/bin/env python
'''
============================================================
Hutt Kids - compute chr 22 phasing stats. 

Created on July 15, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os, numpy as np
from impute.imputation.ImputationSet import ImputationSet

path = os.environ['OBER_OUT'] + '/kids'
chrom = 22
prefix = path + '/cytosnp/chr%d/cytosnp.imputed' % (chrom,)

# Load phasing result
problem = im.io.read_npz(prefix + '.phased.npz')
t = ImputationSet.from_problem(problem)
c = im.v.iv.stats_computer(t, problem.g)
s = c.stats()

