'''
============================================================
Test the IBD posterior obtained by HMM-haplotype in isolation
for several simple cases. 

Created on January 24, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, hmmt, itemutil
from numpy.ma.testutils import assert_almost_equal, assert_equal
from impute.ibd.distant.ibd_hmm_hap import ProbIbdHmmHapCalculator

class TestIbdHmmHap(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        np.set_printoptions(threshold=np.nan, linewidth=200)

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_symmetry_small_lam(self):
        '''Test the HMM on simulated sib data in a fully symmetric setting. The posterior must also
        be symmetric (note: the Delta used here is unrealistic).'''
        self._test_symmetry(0.5, 0.25, 10)

    def test_symmetry_large_lam(self):
        self._test_symmetry(1.0, 0.25, 10)

    def test_symmetry_large_lam_small_m(self):
        self._test_symmetry(1.0, 0.25, 2, e=0)

    def test_likelihood_small_m(self):
        '''Test liklihood of symmetric and non-symmetric paths for 2 SNPs.'''
        m = 2
        c = self._hmm_calculator(1.0, 0.25, m, e=0)
        Obs = np.tile(3, (m,))
        Q = list(itemutil.binseq(m))
        Q_star = Q[np.argmax(np.array([hmmt.log_likelihood(c.m, Obs, x) for x in Q]))]
        c.prob()
        assert_equal(c.Q_star, Q_star, 'Wrong Viterbi path (when compared with Brute-force maximization)')
        
            
    def test_sibs(self):
        '''Test the HMM on simulated sib data.'''
        m = 10  # no. of SNPs
        prob = ProbIbdHmmHapCalculator(lam=0.5, f=0.25,
                                       x=np.arange(m),
                                       p=0.5 * np.ones(m),
                                       h1=np.ones((m,), dtype=np.uint),
                                       h2=np.ones((m,), dtype=np.uint),
                                       e=0.01).prob()
        # print prob.tolist()
        assert_almost_equal(prob,
                            [0.6608271143031219, 0.7480821152397511, 0.7945731630815513, 0.8182183564299464, 0.828130781040886, 0.8281307810408859, 0.8182183564299463, 0.7945731630815512, 0.7480821152397511, 0.660827114303122],
                            decimal=5, err_msg='Wrong posterior IBD')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def _hmm_calculator(self, lam, f, num_snps, e=0.01):
        '''Return the HMM calculator for this test case.'''
        m = num_snps
        return ProbIbdHmmHapCalculator(lam=lam, f=f,
                                       x=np.arange(m),
                                       p=0.5 * np.ones(m),
                                       h1=np.ones((m,), dtype=np.uint),
                                       h2=np.ones((m,), dtype=np.uint),
                                       e=e)
    
    def _test_symmetry(self, lam, f, num_snps, e=0.01):
        '''Test the HMM on simulated sib data in a fully symmetric setting. The posterior must also
        be symmetric (note: the Delta used here is unrealistic).'''
        calculator = self._hmm_calculator(lam, f, num_snps, e)
        prob = calculator.prob()
        Q_star = calculator.Q_star
        assert_almost_equal(prob, np.flipud(prob),
                            decimal=5, err_msg='Posterior IBD not symmetric in a symmetric setting')
        assert_equal(Q_star, np.flipud(Q_star),
                     err_msg='Viterbi path not symmetric in a symmetric setting')
