#!/usr/bin/env python
'''
============================================================
Imputation test - chromosome 22.

Created on February 4, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, os

OBER = os.environ['OBER']
p = im.hutt('hutt.phased.npz')
ibd = im.smart_segment_set.SmartSegmentSet.load(p.pedigree.num_genotyped, OBER + '/out/segments.out')
t = im.imputation.ImputationSet.from_file(p.pedigree, OBER + '/data/impute/rare/rare.npz') #@UndefinedVariable
snps = np.where(t.snp['chrom'] == 22)[0] # SNP list out of all SNPs in t to impute 
im.imputation.iibd.impute(p.haplotype, ibd, t, snp=snps, debug=False)
