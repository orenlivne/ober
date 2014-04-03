'''
============================================================
Test imputation of a single Affy SNP on chr22 (validation).

Created on March 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest#, impute as im, numpy as np
#from numpy.ma.testutils import assert_equal, assert_almost_equal
#from impute.impute_test_util import assert_size_equals

class TestImputationSnp(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        # The way to load a pedigree in conjunction with a genotype set is to recode
        # its sample IDs to consecutive for easier access by phasers.

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_snp_1200(self):
        '''Test imputation accuracy and call rate at SNP 1200, a typical SNP in the middle of chr22
        where the IBD segments probably have ideal coverage.'''
        pass

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
