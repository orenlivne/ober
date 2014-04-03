'''
============================================================
Test overriding imputation results with a plink file of
Affymetrix/IMPUTE2 data.

Created on December 17, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, impute as im
from impute.cgi.override_imputed_by_plink import override_imputed_by_plink
from numpy.ma.testutils import assert_equal

class TestOverrideImputedByPlink(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        # Test file directory
        self.dir = im.itu.OVERRIDE_IMPUTED_BY_PLINK

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_override_with_conordant_data(self):
        '''Test overriding an imputed result by a concordant PLINK data.'''
        
        overriden_line = list(override_imputed_by_plink(self.dir + '/chr7_mutation.tsv',
                                                        self.dir + '/chr7_mutation.impute2',
                                                        override_tag=im.constants.LD_IMPUTED))[0]
        expected = np.loadtxt(open(self.dir + '/chr7_mutation.impute2.out', 'rb'), dtype=str)
        assert_equal(overriden_line, expected, 'Wrong overridden genotypes')

    def test_override_with_discordant_data(self):
        '''Test overriding an imputed result by a discordant PLINK data. Should not override genotypes.'''
        overriden_line = list(override_imputed_by_plink(self.dir + '/chr7_mutation.tsv',
                                                        self.dir + '/chr7_mutation.impute2.discordant',
                                                        override_tag=im.constants.LD_IMPUTED,
                                                        max_mismatch_partial=0))[0]
        expected = np.loadtxt(open(self.dir + '/chr7_mutation.impute2.discordant.out', 'rb'), dtype=str)
        assert_equal(overriden_line, expected, 'Wrong overridden genotypes')
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
