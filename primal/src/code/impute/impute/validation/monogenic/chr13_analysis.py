#!/usr/bin/env python
'''
============================================================
Candadian Hutterites SNP problem: 8 discordant individuals
due to the IMPUTE2 LD-based method are imputed 11 instead of
00. But we input our haps into it and it seems that they are
IBS with some WGS carriers around the mutation. 

Created on December 16, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os, numpy as np, itertools as it  # , matplotlib.pyplot as P

def are_discordant_ibs_with_wgs_carriers(data_dir):
    chrom, bp = 13, 73539488
    
    p = im.io.read_npz('%s/phasing/chr%d/hutt.phased.npz' % (os.environ['OBER_OUT'], chrom))
    # W = im.wgs_sample_index()  # @UndefinedVariable
    W_carrier = [p.pedigree.node_of[x] for x in np.loadtxt(data_dir + '/wgs-carriers', dtype=int)]  # @UndefinedVariable
    I = [p.pedigree.node_of[x] for x in np.loadtxt(data_dir + '/imputed-carriers', dtype=int)]  # @UndefinedVariable
    
    # Compare discorant individuals' haplotypes with WGS carriers' aroudn the mutation
    s = p.haplotype.nearest_snp(bp)
    t = 11
    snps = np.arange(s - t, s + t + 1)
    
    for i in it.islice(I, 0, 1):
        print 'Discordant individual %d [%d]: matching WGS carrier haplotypes' % (p.pedigree.sample_id[i], i)
        for x in ((w, p.pedigree.sample_id[w], d) for w in W_carrier 
                  for d in im.diff.all_diffs(p.h[snps], i, w) 
                  if sum(abs(d)) <= 3):
            print '%-8d [%-4d] %s' % (x[1], x[0], ' '.join(map(str, x[2])))

####################################################################################
if __name__ == '__main__':
    data_dir = os.environ['OBER_OUT'] + '/requests/taqman'

    # are_discordant_ibs_with_wgs_carriers(data_dir)
    # problem = im.io.read_plink(prefix=data_dir + '/taqman', pedigree=im.itu.HUTT_PED, haplotype=None, frames=None)
    # p, t = im.v.iv.impute_problem(problem, debug=2, debug_sample=1379)
    
    problem = im.io.read_plink(prefix=data_dir + '/taqman.genotypes', pedigree=im.itu.HUTT_PED, haplotype=None, frames=None, pedigree_genotyped=im.itu.HUTT_GENOTYPED_PED)
    p, t = im.v.iv.impute_problem(problem, index_file=data_dir + '/taqman.genotypes.index', genotypes_as_is=True, debug=2)
