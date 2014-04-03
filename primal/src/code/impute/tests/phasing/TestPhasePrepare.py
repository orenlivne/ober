'''
============================================================
Test phasing pre-processing.

Created on October 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np
from impute import impute_test_util as itu
from impute.data import io
from impute.phasing.pre_processing import prepare_phaser
from numpy.ma.testutils import assert_almost_equal

class TestPhasePrepare(unittest.TestCase):
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
        self.phaser = prepare_phaser()

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_estimate_genotype_frequencies(self):
        '''Test estimating genotype frequencies for each SNP using the prepare module.'''
        problem = io.read_plink(pedigree=itu.HUTT_PED, prefix=itu.GENOTYPE_SAMPLE, haplotype=None)
        itu.assert_size_equals(problem.genotype, 8, 1415)
        self.phaser.run(problem)
        # Expected result
        frequency = np.array([[  3.44876319e-01,   4.74911660e-01,   1.75265014e-01,   4.94699646e-03],
                              [  1.06007066e-02,   1.86572433e-01,   8.02826881e-01,   0.00000000e+00],
                              [  1.83745585e-02,   2.93286204e-01,   6.63604259e-01,   2.47349832e-02],
                              [  2.17667848e-01,   5.17314494e-01,   2.64310956e-01,   7.06713763e-04],
                              [  2.16961130e-01,   5.16607761e-01,   2.56537110e-01,   9.89399292e-03],
                              [  6.10600710e-01,   3.33568901e-01,   5.58303893e-02,   0.00000000e+00],
                              [  7.24381626e-01,   1.59010604e-01,   7.06713763e-04,   1.15901060e-01],
                              [  4.19081271e-01,   4.62190807e-01,   1.11660779e-01,   7.06713786e-03]])
        assert_almost_equal(problem.info.genotype_frequency, frequency,
                            decimal=5, err_msg='Wrong SNP genotype frequency estimation')
#        assert_almost_equal(problem.info.allele_frequency(1), frequency[:, 0] + 0.5 * frequency[:, 1],
#                            decimal=5, err_msg='Wrong SNP genotype allele frequency estimation')
#        assert_almost_equal(problem.info.allele_frequency(2), frequency[:, 2] + 0.5 * frequency[:, 1],
#                            decimal=5, err_msg='Wrong SNP genotype allele frequency estimation')
        assert_almost_equal(problem.info.allele_frequency(2), 1.0 - problem.info.allele_frequency(1),
                            decimal=5, err_msg='Wrong SNP genotype allele frequency estimation')
        #itu.assert_problem_stats(problem, 22640, 18567, 10)
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
