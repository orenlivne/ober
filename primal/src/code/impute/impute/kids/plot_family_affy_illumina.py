#!/usr/bin/env python
'''
============================================================
Estimate Hutterite kid imputation rates so that Carole can
decide whether to genotype them with a dense or sparse
Illumina chip.   

Created on July 15, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os, numpy as np, matplotlib.pyplot as P

# Load data
ped = im.hutt_pedigree()
path = os.environ['OBER_OUT'] + '/kids'
chrom = 22
prefix = path + '/cytosnp/chr%d/cytosnp.imputed' % (chrom,)
illumina = im.io.read_npz(prefix + '.phased.npz')
affy = im.hutt('hutt.phased.npz')

# Large family - with lots of sibs of one of the new Hutt kids
#parents = 246, 389
parents = 288, 465
f = ped.find_family(parents[0], parents[1])

# Compare Illumina, Affy IBD sharing pictures
P.figure(1)
im.plots.plot_family_comparison(affy, f, 1, xaxis='bp')
P.savefig(os.environ['OBER'] + '/doc/kids/family_%d_%d_affy.png' % parents)

P.figure(2)
im.plots.plot_family_comparison(illumina, f, 1, xaxis='bp')
P.savefig(os.environ['OBER'] + '/doc/kids/family_%d_%d_illumina.png' % parents)
