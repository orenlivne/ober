#!/usr/bin/env python
'''
============================================================
Feasibility study of imputing untyped ancestors using
POO information on quasi-founders.
determination.

Created on August 15, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, util, matplotlib.pylab as P, os, numpy as np, itertools as it
#@from db_gene.snp.file_dao import ChromDao
from impute.data.constants import MATERNAL, PATERNAL#, CHROMOSOMES
from collections import Counter

def analyze_family(p, samples, min_length=0, max_colors=0, title='', plot=True, debug=False):
    '''Analyze QF poo-aligned sibs, get IBD segments of all paternal & maternal haps.''' 
    # segments = im.segment.SegmentSet.load(open('%s/poo/segments-family5.out' % (os.environ['OBER_OUT'],), 'rb'))
    segments = im.ibd_distant_hap.among_samples_segments(p, samples, im.phase.PhaseParam())
    segments = segments.segments_larger_than(min_length)
    pa = im.color.hap_color.hap_colors(list(it.product(sorted(samples), im.constants.ALLELES)), segments, max_colors=4, max_sweeps=3)
    if plot: im.color.hap_color.plot_palette(pa, pair_gap=10, linewidth=0.2, title=title, show_primary_haps=True, show_regions=debug)
    return pa, segments

def family_haps(chrom, ibd, s):
    '''Retrieve the <=four haplotypes of a set of sibs from the IBD segment index.'''
    region_size = int(ibd._region_size)
    ibd._load_chrom(chrom)  # Force loading of this chromosome
    num_snps = len(ibd._snp)
    group = np.zeros((num_snps, len(s), 2), dtype=np.int16)
    print group.shape
    print ibd
    for start in xrange(0, num_snps, region_size):
        ibd.find(chrom, start, 0, 0)  # Force loading of this region
        region_data = ibd._group_index[:, s, :]
        print ibd._snp
        print ibd._chrom
        print region_data.shape
        print group.shape
        group[start:start + region_data.shape[0]] = region_data
    return group

def plot_hist_num_sibs(num_sibs):
    '''Histogram of #QF sibs.'''
    P.figure(1)
    P.clf()
    P.hist(num_sibs)    
    P.title('# Typed Sibs of Quasi Founders')
    P.ylabel('Frequency')
    P.xlabel('# Sibs')
    # P.legend(loc='upper right', prop={'size': 10})
    P.show()

def parent_coverage_fraction(pa, p):
    '''Return the % of base pairs out of the base pairs covered by p\'s SNPs for which we
    can reconstruct each of the four parental haplotypes.''' 
    return np.array(map(pa.parent_coverage_bp, xrange(4))) / float(p.info.snp['base_pair'][-1] - p.info.snp['base_pair'][0])

def qf_families_parent_coverage(p):
    '''Analyze all QF families in the Problem object p. Return (parents, parent hap bp coverage) for all families.'''
    d = Counter()
    for f in frozenset(p.pedigree.find_family_by_child(i, genotyped=False) for i in p.pedigree.quasi_founders):
        children = im.gt.genotyped_children(p, f)
        if len(children) > 1:
            pa = analyze_family(p, children, plot=False)[0]
            a, b = len(children), np.array(map(pa.parent_coverage_bp, xrange(4)))  # parent_coverage_fraction(pa, p)
            _, _, paternal_colors, maternal_colors = im.color.hap_color.best_hap_alignment_to_colors(pa)
            print f, a, parent_coverage_fraction(pa, p), paternal_colors, maternal_colors
            parental_colors = (paternal_colors, maternal_colors)
            for parent_type in (PATERNAL, MATERNAL):
                for i in xrange(2): d[(f.parents[parent_type], i)] = b[parental_colors[parent_type][i]]
    return d

def plot_two_families():
    '''Test ancestor imputation and child POO alignment for two families on chromosome 22.'''
    # Parameters
    chrom = 22
    plot = False  # True
    save_plot = False  # True
    debug = False  # True
    
    # Read data
    p = im.hutt('hutt.phased.npz')
    q = p.pedigree.quasi_founders
    # aligned = set(p.haplotype.aligned_samples)
    t = frozenset([frozenset(im.gt.genotyped_children(p, p.pedigree.find_family_by_child(i, genotyped=False))) for i in q])
    num_sibs = map(len, t)
    print 'Distribution of QF family sizes', util.occur_dict(num_sibs)
    # plot_hist_num_sibs(num_sibs)
    
    # ibd = im.index.segment_index.SegmentIndex(os.environ['OBER_OUT'] + '/index_segments')
    
    if plot: P.figure(1)
    s = set([x for x in t if 1049 in x][0]) - set([1049])
    pa, _ = analyze_family(p, s, max_colors=4, title='Haplotype Coloring: Quasi-Founder Sibs, All, Chrom. %d' % (chrom,), plot=plot, debug=debug)
    if save_plot: P.savefig(os.environ['OBER'] + '/doc/poo/qf_family/hap_colors_poo.png')
    
    if plot: P.figure(2)
    s2 = set([x for x in t if 1049 in x][0])
    analyze_family(p, s2, max_colors=4, title='Haplotype Coloring: Quasi-Founder Sibs, POOC Chrom. %d' % (chrom,), plot=plot, debug=debug)
    if save_plot: P.savefig(os.environ['OBER'] + '/doc/poo/qf_family/hap_colors_all.png')
    
    if plot: P.figure(3)
    f = p.find_family(10, 1414)  # 4 children, genotyped parents
    s3 = f.children
    analyze_family(p, s3, max_colors=4, title='Haplotype Coloring: Non-Founder Sibs Chrom. %d' % (chrom,), plot=plot, debug=debug)
    if save_plot: P.savefig(os.environ['OBER'] + '/doc/poo/nf_family/hap_colors.png')
    
    if plot: P.show()
    print im.color.hap_color.best_hap_alignment_to_colors(pa)
    print 'Regions', pa.num_regions
    print 'Parental haplotype coverage %', parent_coverage_fraction(pa, p)
    print 'Children coverage by parental haplotypes', pa.color_sequence_coverage(np.arange(4))

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
#     parent_coverage = Counter()
#     for chrom in CHROMOSOMES[-1:]:
#         print 'Chromosome', chrom
#         p = im.io.read_npz('%s/phasing/chr%d/hutt.phased.npz' % (os.environ['OBER_OUT'], chrom))
#         # plot_two_families()
#         parent_coverage_chrom = qf_families_parent_coverage(p)
#         parent_coverage += parent_coverage_chrom 
#         print parent_coverage_chrom
#     parents = set(x[0] for x in parent_coverage.iterkeys())
#     a = [(b, (parent_coverage[(b, 0)] + parent_coverage[(b, 1)]) / (2.*sum(ChromDao.TOTAL_BP_TYPED[-1:]))) for b in set(x[0] for x in parent_coverage.iterkeys())]
#    print a

    p = im.hutt('hutt.phased.npz')

    # pa, segments = analyze_family(p, np.setdiff1d(im.gt.genotyped_children(p, f), [1069]))#, debug=True)
    # f = p.find_family_by_child(998, genotyped=False)
    # pa, segments = analyze_family(p, im.gt.genotyped_children(p, f))
    
    f = p.find_family_by_child(640, genotyped=False)
    pa, segments = analyze_family(p, im.gt.genotyped_children(p, f))
