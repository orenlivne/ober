#!/usr/bin/env python
'''
============================================================
Run stage 3 on families around sample 507. 

Created on September 13, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, operator, itertools
from impute import impute_test_util as itu
from impute.ibd import ibd_germline as ig
#import numpy as np
from impute.data.constants import ALLELES

def compress(segments, haps):
    d = {}
    for s in segments: 
        d.setdefault(s.snp, []).append(s.samples)
    for ibd in d.itervalues():
        for remaining in haps - reduce(set.union, ibd):
            ibd.append(set([remaining]))
    return d

def longest_segment(d): 
    return max([(snp, len(ibd)) for (snp,ibd) in d.iteritems()], key=operator.itemgetter(1))

def hash_ibd_list(ibd):
    return dict(((sample, i) for (i, samples) in enumerate(ibd) for sample in samples))

def highest_degree_node(ibd):
    return max([(i,len(ibd[hap_index[(i,0)]])+len(ibd[hap_index[(i,1)]])) for i in sibs], 
               key=operator.itemgetter(1))[0]

p = im.io.read_npz(itu.SIB_FOUNDERS_STAGE3)
f = p.families(genotyped=False)[0]
sibs = ig._filled_members(p, f) 
h_mat = ig._HapMatrix(p, sibs)
print h_mat
c = ig.GermlineIbdComputer()
segments = c.ibd_segments(h_mat)
segments.group_to_disjoint()
print segments

haps = set(itertools.product(sibs, ALLELES))
d = compress(segments, haps)
(best_segment, num_haps) = longest_segment(d)
print best_segment, num_haps
ibd = d[best_segment]
print ibd
hap_index = hash_ibd_list(ibd)
template = highest_degree_node(ibd)
print template

#stats = np.array([(len(s.samples), s.length) for s in m])
#best_segment = np.argsort(np.lexsort((-stats[:,1],-stats[:,0])))[0]










