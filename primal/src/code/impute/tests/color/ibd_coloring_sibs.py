#!/usr/bin/env python
'''
============================================================
Test IBD plotting in a nuclear family. 

Created on January 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, itertools, matplotlib.pyplot as P, sys
from impute.tools.param import PhaseParam

# print sys.argv
generate_plots = True if (len(sys.argv) < 2) else bool(int(sys.argv[1]))
p = im.io.read_npz(im.itu.FAMILY_TOO_ZEROED_STAGE2) 
haps = list(itertools.product(im.gt.genotyped_members(p, p.first_family), xrange(2)))
children = im.gt.genotyped_children(p, p.first_family)

# IBD between sib pairs
child = 3
sib = 2
sib_ibd = im.ibd_distant.ibd_segments_with_relatives(p, child, [sib],
                                                     PhaseParam(margin=0., surrogate_parent_fill_threshold=0.9, debug=True),
                                                     im.ibd_hmm.prob_ibd_hmm)
print sib_ibd
if generate_plots:
    P.figure(1)
    child_haps = list(itertools.product([child, sib], xrange(2)))
    g = im.plots.plot_hap_coloring(sib_ibd, child_haps, pair_gap=10, linewidth=6, title='Sib IBD Segments',
                                   snp_range=(0, p.num_snps))
#    P.savefig(os.environ['OBER'] + '/doc/ibd/hmm/family_ibd_hmm_sib.png')
