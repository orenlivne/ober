#!/usr/bin/env python
'''
============================================================
Test GERMLINE IBD on 507's ungenotyped family. 

Created on September 15, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im
import numpy as np

p = im.hutt('hutt.stage3.npz')
q = im.hutt('hutt.stage3.npz')
phaser = im.phase_distant.family_sib_comparison_phaser()
i = 507
phaser.run(q, im.PhaseParam(single_member=i, debug=True))

print np.where(p.haplotype.data[:,i,:] != q.haplotype.data[:,i,:])
