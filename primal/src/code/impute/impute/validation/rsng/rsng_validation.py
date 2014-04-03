'''
Created on May 3, 2013

@author: oren
'''
import impute as im, os

# problem = im.io.read_plink(prefix=os.environ['OBER_OUT'] + '/requests/rsng/plink/rsng', pedigree=im.itu.HUTT_PED, haplotype=None, frames=None)
# im.io.write_npz(problem, os.environ['OBER_OUT'] + '/requests/rsng/plink/rsng.npz')
problem = im.io.read_npz(os.environ['OBER_OUT'] + '/requests/rsng/plink/rsng.npz')
ibd_sample_index = im.io_pedigree.genotyped_sample_index(im.examples.hutt_pedigree(), problem.pedigree)
p, t = im.v.iv.impute_problem(problem, snps=[46], debug=2, ibd_sample_index=ibd_sample_index)
