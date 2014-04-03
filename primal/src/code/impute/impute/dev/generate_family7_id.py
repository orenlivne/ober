'''
============================================================
Generate the family 7 identity coefficient file.

Created on January 21, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, db_gene, itertools

# Temporarily comment out the frames field in im.io functions for this cmomand to work
p = im.io.read_npz(im.itu.FAMILY7 + '.npz')

i = p.pedigree.genotyped_sample_id()
k = db_gene.snp.file_dao.IdCoefDao(im.PhaseParam().id_coef_file)
a = np.array([(x, y, m[0]) + tuple(m[1].tolist()) for (x, y, m) in  ((x, y, k.id_coefs(x, y)) for x, y in itertools.product(i, i))])
np.savetxt(im.itu.FAMILY7 + '.id', a, fmt='%d %d %e %e %e %e %e %e %e %e %e %e')
