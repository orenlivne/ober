#!/usr/bin/env python
'''
============================================================
Subjects 80, 248 are sibs in the pedigree but might be only
half-sibs. Same mom different dads.
Find whether 80 and its dad (547) share IBS>=1 everywhere.

Created on December 6, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
from __future__ import division
import impute as im, os , numpy as np
from impute.data.constants import CHROMOSOMES

####################################################################################
if __name__ == '__main__':

    proband, parent = 579,1403 #248, 547#80, 547
    
    yes_ibs, no_ibs = 0, 0
    for chrom in CHROMOSOMES:
        p = im.io.read_npz('%s/phasing/chr%d/hutt.npz' % (os.environ['OBER_OUT'], chrom))
        ibs = im.diff.ibs_state(p.g, proband, parent)
        no, yes = len(np.where(ibs == 0)[0]), len(np.where(ibs > 0)[0])
        print 'Chromosome %d, %%IBS>=1 = %f' % (chrom, yes / (yes + no))
        no_ibs += no
        yes_ibs += yes
    print 'Genome %%IBS>=1 = %f' % (yes_ibs / (yes_ibs + no_ibs),)
