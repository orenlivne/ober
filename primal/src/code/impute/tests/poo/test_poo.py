#!/usr/bin/env python
'''
============================================================
Test IBD clique kinship comparison for parent-of-origin
determination.

Created on July 2, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, time, numpy as np, os

'''
--------------------------------------------------
Main program
--------------------------------------------------
'''
chrom = 22
# a = im.poo.Aligner(chrom)

# print a.debug_sample_snp(1101, 1500)
# print a.debug_sample(6)
# a.debug_sample_snp(6, 2000)
# a.debug_sample_snp(0, 500)

# Time POO vs. # processes used
for num_processes in [1]:#[2, 4]:
    start_time = time.time()
    m = im.poo.determine_poo(chrom, params=im.PhaseParam(debug=True, poo_snp_step_size=100, num_processes=num_processes))
    print '#processes=%d, time %.2f' % (num_processes, time.time() - start_time)
 
np.savetxt(os.environ['OBER'] + '/doc/poo/m-chr%d.txt' % (chrom,), m)
im.poo.plot_flip_measure_vs_sample(chrom, m)
