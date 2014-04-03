#!/usr/bin/env python
'''
============================================================
Test IBD plotting in a nuclear family. 

Created on January 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, itertools, matplotlib.pyplot as P, sys
from impute.color.hap_color_grouping import plot_hap_coloring

print sys.argv
generate_plots = True if (len(sys.argv) < 2) else bool(int(sys.argv[1]))
p = im.io.read_npz(im.itu.FAMILY_TOO_ZEROED_STAGE4) 
haps = list(itertools.product(im.gt.genotyped_members(p, p.first_family), xrange(2)))
children = im.gt.genotyped_children(p, p.first_family)
child_haps = list(itertools.product(children, xrange(2)))

params = im.PhaseParam()#debug=True)

#-------------------------------------------------------
# Nuclear family phasing
#-------------------------------------------------------
ibd = p.info.ibd
print ibd
if generate_plots:
    # Including parents
    P.figure(1)
    plot_hap_coloring(ibd, haps, pair_gap=10, linewidth=6, title='Family IBD Segments: Nuclear Family Phasing')
    #P.savefig(os.environ['OBER'] + '/doc/ibd/hmm/family_ibd_nuclear.png')
    # Nuclear family phasing without parents
#    P.figure(2)
#    plot_hap_coloring(ibd, child_haps, pair_gap=10, linewidth=6, title='Family IBD Segments: Nuclear Family Phasing')
    #P.savefig(os.environ['OBER'] + '/doc/ibd/hmm/family_ibd_nuclear.png')

#-------------------------------------------------------
# GG-IBD between each parent and all children 
#-------------------------------------------------------
distant_ibd = sum((im.ibd_distant.ibd_segments_with_relatives(p, parent, children,
                                                              params,
                                                              im.ibd_hmm.prob_ibd_hmm)
                   for parent in p.first_family.parents), im.segment.SegmentSet([]))
print distant_ibd
if generate_plots:
    P.figure(3)
    g = plot_hap_coloring(distant_ibd, haps, pair_gap=10, linewidth=6, title='Family IBD Segments: HMM, Parent vs. Children')
    #P.savefig(os.environ['OBER'] + '/doc/ibd/hmm/family_ibd_hmm.png')

#-------------------------------------------------------
# GG-IBD between sib pairs
#-------------------------------------------------------
sib_ibd = sum((im.ibd_distant.ibd_segments_with_relatives(p, child, list(set(children) - set([child])),
                                                          params,
                                                          im.ibd_hmm.prob_ibd_hmm)
               for child in children), im.segment.SegmentSet([]))
print sib_ibd
if generate_plots:
    P.figure(4)
    g = plot_hap_coloring(sib_ibd, child_haps, pair_gap=10, linewidth=6,
                                   title='Family IBD Segments: HMM, All Sib Pairs', snp_range=(0, p.num_snps))
    #P.savefig(os.environ['OBER'] + '/doc/ibd/hmm/family_ibd_hmm_sib.png')

#-------------------------------------------------------
# HH-IBD between sib pairs
#-------------------------------------------------------
sib_ibd_hh = im.ih.among_samples_segments(p, children, params)
print sib_ibd_hh
if generate_plots:
    P.figure(5)
    g = plot_hap_coloring(sib_ibd_hh, child_haps, pair_gap=10, linewidth=6,
                                   title='Family IBD Segments: HMM-HH, All Sib Pairs', snp_range=(0, p.num_snps))
    #P.savefig(os.environ['OBER'] + '/doc/ibd/hmm/family_ibd_hmm_sib.png')
