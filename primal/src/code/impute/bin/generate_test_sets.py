#!/usr/bin/env python
'''
============================================================
Re-generate all test sets (chr22 data).

Created on August 16, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
from impute.phasing import examples
from impute.preprocess import convert
from impute.phasing.examples import CHR22
from impute import impute_test_util as itu

# Need to be done once given the original PLINK TPED Hutterites data set
convert.main(pedigree=itu.HUTT_PED, prefix=CHR22 + '/hutt', npz=CHR22 + '/hutt.npz', target='npz', debug=True)

for stage in xrange(1, 5):
    examples.phase_chr22('hutt.%snpz' % ('stage%d.' % (stage - 1,) if stage > 1 else ''),
                         'hutt.stage%d.npz' % (stage,),
                         stage=stage)

examples.generate_test_sets()
