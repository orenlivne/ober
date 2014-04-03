#!/usr/bin/env python
'''
============================================================
Distant IBD (using IBS for now) for the sib founder family
of 507. 

Created on December 10, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im  # , numpy as np
from impute import impute_test_util as itu, ibd
from impute.ibd.ibd_distant import IbdDistantSegmentCalculator
# from impute.data.constants import MEGA_BASE_PAIR
# from impute.ibd.ibd import _phase_in_ibd_segment
# from impute.tools.param import PhaseParam

# def compress(segments, haps):
#    d = {}
#    for s in segments: 
#        d.setdefault(s.snp, []).append(s.samples)
#    for ibd in d.itervalues():
#        for remaining in haps - reduce(set.union, ibd):
#            ibd.append(set([remaining]))
#    return d
#
# def longest_segment(d): 
#    return max([(snp, len(ibd)) for (snp,ibd) in d.iteritems()], key=operator.itemgetter(1))
#
# def hash_ibd_list(ibd):
#    return dict(((sample, i) for (i, samples) in enumerate(ibd) for sample in samples))
#
# def highest_degree_node(ibd):
#    return max([(i,len(ibd[hap_index[(i,0)]])+len(ibd[hap_index[(i,1)]])) for i in sibs], 
#               key=operator.itemgetter(1))[0]

p = im.io.read_npz(itu.SIB_FOUNDERS_STAGE3)
print p.fill_fraction()
f = p.families(genotyped=False)[0]
i = 2
print p.fill_fraction()

calculator = IbdDistantSegmentCalculator(p, max_path_length=2, debug=False)
#fill_threshold=0.9 