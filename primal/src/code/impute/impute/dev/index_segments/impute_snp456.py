#!/usr/bin/env python
'''
============================================================
Impute SNP 456 on chr22.

Created on April 1, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os

def run_snp(snp, min_len=0.0):
    im.index.index_segments.main(None, im.examples.CHR22 + '/hutt.phased.info.npz', im.examples.CHR22 + '/segments.out', os.environ['OBER_OUT'] + '/ibd/chr22/index_segments_%d' % (snp,), num_processes=1, region_size=1, snp_index=snp, min_len=min_len, debug=2, algorithm='amg')
    return im.v.iv.impute_chrom(22, snps=[snp], debug=2, input_dir=os.environ['OBER_OUT'] + '/phasing/chr22', segment_file=os.environ['OBER_OUT'] + '/ibd/chr22/index_segments_%d' % (snp,), algorithm='index')

p, t, i = run_snp(3000)
