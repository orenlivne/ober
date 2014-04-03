#!/usr/bin/env python
'''
============================================================
Extract family IBD sharing picture for an individual which
Rebecca inquired about (at the monogenic chr11 mutation).

Created on July 22, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os, matplotlib.pyplot as P, numpy as np

data_dir = os.environ['OBER_OUT'] + '/requests/monogenic/chr11_88911393'
out_dir = os.environ['OBER'] + '/doc/imputation/monogenic'

# Extract data on family of suspected carrier individual 512. Father: 592, mother: 46
problem = im.io.read_npz(os.environ['OBER_OUT'] + '/phasing/chr11/hutt.phased.npz')
f = problem.find_family(592, 46)
ibd = im.index.segment_index.SegmentIndex('/home/oren/ober/out/index_segments/')
snp = ibd.nearest_left_snp(11, 88911393)
print 'snp', snp

# IBD picture of 512's nuclear family
d = im.plots.plot_family_comparison(problem, f, 0, yaxis='id')
P.savefig(os.environ['OBER'] + '/doc/imputation/monogenic/chr11_family_paternal_ibd.png')
# Draw a vertical line at the location of the mutation
P.plot([snp, snp], [-1, d[0].shape[1]], 'k-')

# Validate that IBD segments are consistent with that picture at the Affy snp closest (on the left) to the mutation 
for allele in im.constants.ALLELES:
    print 'Children haps flagged in IBD cliques as sharing the father''s %s haplotype:' % (im.constants.ALLELE_LABELS[allele])
    print [tuple(x) for x in ibd.find(11, 8602, 592, allele) if x[0] in f.children]

# Draw pedigree of imputed genotypes of father's family, where we think the 'A' allele was inherited to 512 from
g = np.loadtxt(data_dir + '/chr11_88911393.lgen', usecols=[1, 2, 3], dtype=[('id', 'i4'), ('a1', 'S12'), ('a2', 'S12')])
labels = dict((i, a1 + '|' + a2) for i, a1, a2 in g)
i = 592

#nbhrs = problem.pedigree.neighbors_genotyped_selective(i, 2)
#nbhrs.append(1489)
#nbhrs = np.array(list(set([x for f in set([f for i in (x for x in nbhrs if not set(problem.pedigree.graph.predecessors(x)) & set(nbhrs)) for f in problem.pedigree.find_families_by_father(i, genotyped=False) + problem.pedigree.find_families_by_mother(i, genotyped=False)]) for x in f.member_list])))
#pedigree = problem.pedigree.sub_pedigree(nbhrs)
#im.pt.to_pedfiddler(pedigree, data_dir + '/father_family.dat', labels=labels)
#im.pt.draw_pedigree(pedigree, out_dir + '/chr11_father_pedigree.png', labels=labels)
p = im.pt.draw_member_neighbor_genotyped_pedigree(problem, i, 2, out_dir + '/chr11_father_pedigree.png', labels=labels)

# IBD picture of 592's nuclear family
f = problem.find_family(1753, 944, genotyped=False)
d = im.plots.plot_family_comparison(problem, f, 1, yaxis='id')
# Draw a vertical line at the location of the mutation
P.plot([snp, snp], [-1, d[0].shape[1]], 'k-')
P.savefig(out_dir + '/chr11_father_family_paternal_ibd.png')

s = im.examples.problem_ibd_segments(problem, 944, None, 592, 1, debug=False)
