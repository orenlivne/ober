'''
============================================================
Test the IBD posterior obtained by HMM in isolation for
several simple cases. 

Created on January 7, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np
from impute.ibd.distant.ibd_hmm import ProbIbdHmmCalculator
from numpy.ma.testutils import assert_almost_equal, assert_equal

class TestIbdHmm(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    # Identity coefficients between outbred sibs
    DELTA_SIBS = np.array([0, 0, 0, 0, 0, 0, 0.25, 0.5, 0.25])
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_symmetry_small_lam(self):
        '''Test the HMM on simulated sib data in a fully symmetric setting. The posterior must also
        be symmetric (note: the Delta used here is unrealistic).'''
        self._test_symmetry(0.5, 10)

    def test_symmetry_large_lam(self):
        self._test_symmetry(1.0, 2, e=0)

    def test_sibs(self):
        '''Test the HMM on simulated sib data.'''
        m = 10  # no. of SNPs
        prob = ProbIbdHmmCalculator(lam=1.0, Delta=TestIbdHmm.DELTA_SIBS,
                                    x=np.arange(m),
                                    p=0.5 * np.ones(m),
                                    g1=np.ones((m, 2), dtype=np.uint),
                                    g2=np.ones((m, 2), dtype=np.uint),
                                    e=0.01).prob()
        assert_almost_equal(prob,
                            [0.9250266434221001, 0.9469142114854412, 0.9519559084940241, 0.9533403567284999, 0.9537394679300106, 0.9537394679300106, 0.9533403567284999, 0.9519559084940241, 0.9469142114854413, 0.9250266434221001],
                            decimal=5, err_msg='Wrong posterior IBD')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def _test_symmetry(self, lam, num_snps, e=0.01):
        '''Test the HMM on simulated sib data in a fully symmetric setting. The posterior must also
        be symmetric (note: the Delta used here is unrealistic).'''
        m = num_snps
        calculator = ProbIbdHmmCalculator(lam=lam, Delta=TestIbdHmm.DELTA_SIBS,
                                          x=np.arange(m),
                                          p=0.5 * np.ones(m),
                                          g1=np.ones((m, 2), dtype=np.uint),
                                          g2=np.ones((m, 2), dtype=np.uint),
                                          e=e)
        prob = calculator.prob()
        Q_star = calculator.Q_star + 1  # convert to 1-based states
        assert_almost_equal(prob, np.flipud(prob),
                            decimal=5, err_msg='Posterior IBD not symmetric in a symmetric setting')
        assert_equal(Q_star, np.flipud(Q_star),
                     err_msg='Viterbi path not symmetric in a symmetric setting')
