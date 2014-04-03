#!/usr/bin/env python
'''
============================================================
Re-generate all Hutt phasing stages' files.

Created on August 16, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im
from impute.preprocess import convert
from impute.phasing.examples import CHR22
from impute.phasing.phase import NUM_STAGES

# Need to be done once given the original PLINK TPED Hutterites data set
p = convert.main(pedigree=im.itu.HUTT_PED, prefix=CHR22 + '/hutt', npz=CHR22 + '/hutt.npz', target='npz', debug=True)

for stage in xrange(1, NUM_STAGES + 1):
    im.examples.phase_chr22('hutt.%snpz' % ('stage%d.' % (stage - 1,) if stage > 1 else ''),
                            'hutt.stage%d.npz' % (stage,), stage=stage)
    # stage=stage, poo_snp_step_size=1600, debug=True, num_processes=2)

im.examples.generate_test_sets()
