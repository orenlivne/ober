#!/usr/bin/env python
'''
============================================================
Debugging family 517 around SNP 2503 on chr5, part 0.
For Carole rare variant imputation grant proposal due
Dec 15.

Stage 4 seems to screw up the family's recombination picture
even though this family should not be modified by this stage\
since both parents are genotyped.

Created on December 7, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, matplotlib.pyplot as P, os

p = im.io.read_npz(os.environ['OBER']+'/out/hutt/phasing.20121207/individual/chr5/hutt.stage3.npz')
f = p.find_family(281, 517)
im.plots.plot_all_family_comparisons(p, f, fig_num=3)
#im.plots.plot_all_family_comparisons(q, f, fig_num=3)
P.show()
