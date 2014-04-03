'''
============================================================
Test imputation statistics functions.

Created on December 5, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, unittest
from numpy.ma.testutils import assert_almost_equal
from impute.imputation import impute_stats_old as istats
from impute.imputation.impute_stats import CALLED_AND_HAS_MINOR_ALLELE, CALLED

class TestImputeStats(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        # Generate two random genotype arrays to compare
        np.random.seed(1)
        (num_snps, num_samples) = (2, 100) 
        self.a = np.sort(np.random.random_integers(0, 2, (num_snps, num_samples, 2)), axis=2) 
        self.b = np.sort(np.random.random_integers(0, 2, (num_snps, num_samples, 2)), axis=2)
        
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_concordance_self(self):
        '''Test concordance between two genotype arrays.'''
        concordance = istats.concordance(self.a, self.a, CALLED)
        TestImputeStats.__assert_concordance_almost_equal(concordance, [1.0, 1.0])

    def test_concordance_with_parameters(self):
        '''Test concordance between two genotype arrays. Apply different parameters.'''
        concordance = istats.concordance(self.a, self.b, CALLED)
        TestImputeStats.__assert_concordance_almost_equal(concordance, [0.45454545, 0.27272727])

        concordance = istats.concordance(self.a, self.b, CALLED, samples=np.arange(90))
        TestImputeStats.__assert_concordance_almost_equal(concordance, [0.5, 0.25])

        concordance = istats.concordance(self.a, self.b, CALLED_AND_HAS_MINOR_ALLELE)
        TestImputeStats.__assert_concordance_almost_equal(concordance, [0.36842105, 0.23809524])
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    @staticmethod
    def __assert_concordance_almost_equal(actual, expected, decimal=5):
        assert_almost_equal(actual[1], expected[1], decimal=decimal, err_msg='Wrong concordance rates')
