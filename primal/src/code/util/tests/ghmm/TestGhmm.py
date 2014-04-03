'''
============================================================
Test the GHMM Hidden Mark Model (HMM) library.

Created on November 20, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, ghmm as g
from UnfairCasino import train_seq 
from numpy.ma.testutils import assert_equal, assert_almost_equal

class TestGhmm(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------

    def setUp(self):
        '''Create a simple dice rolling HMM'''
        self.sigma = g.IntegerRange(1, 7)
        self.A = [[0.9, 0.1], [0.3, 0.7]]
        efair = [1.0 / 6] * 6
        eloaded = [3.0 / 13, 3.0 / 13, 2.0 / 13, 2.0 / 13, 2.0 / 13, 1.0 / 13]
        self.B = [efair, eloaded]
        self.pi = [0.5] * 2
        self.m = g.HMMFromMatrices(self.sigma, g.DiscreteDistribution(self.sigma), self.A, self.B, self.pi)
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_initial_model(self):
        '''Test looking at the initial model; generating observations.'''
        assert_equal(str(self.m),
                     'DiscreteEmissionHMM(N=2, M=6)\n' \
                     '  state 0 (initial=0.50)\n' \
                     '    Emissions: 0.17, 0.17, 0.17, 0.17, 0.17, 0.17\n' \
                     '    Transitions: ->0 (0.90), ->1 (0.10)\n' \
                     '  state 1 (initial=0.50)\n' \
                     '    Emissions: 0.23, 0.23, 0.15, 0.15, 0.15, 0.08\n' \
                     '    Transitions: ->0 (0.30), ->1 (0.70)\n', 'Wrong HMM object representation string')
        TestGhmm.assert_model_matrices_almost_equal(self.m, [self.A, self.B, self.pi])

        # Generate an observation
        size = 20
        obs_seq = self.m.sampleSingle(size)
        obs = map(self.sigma.external, obs_seq)
        assert_equal(len(obs), size, 'Wrong observation vector size')
        
    def test_train_model(self):
        '''Test loading training data and fitting the model parameters.'''
        self.m.baumWelch(train_seq)
        TestGhmm.assert_model_matrices_almost_equal(self.m, 
                                                    [[[0.9042082623992436, 0.09579173760075876], [0.30937033420056315, 0.690629665799436]], 
                                                     [[0.16198458447552336, 0.14425976689332365, 0.17112812060903818, 0.1659152144345754, 0.17723670352872015, 0.17947561005881948], [0.23175180320928723, 0.20230358947748303, 0.16002596788888873, 0.153600271294851, 0.165675856427131, 0.08664251170235841]], 
                                                     [0.48138641074729577, 0.5186135892527041]])
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    @staticmethod
    def assert_model_matrices_almost_equal(m, m_expected, decimal=3):
        '''Check that a model''s matrices equal a set of expected matrices.'''
        (m_A, m_B, m_pi) = m.asMatrices()
        (A, B, pi) = m_expected
        assert_almost_equal(np.array(m_A), np.array(A), decimal=decimal, err_msg='Wrong transition probabilities')
        assert_almost_equal(np.array(m_B), np.array(B), decimal=decimal, err_msg='Wrong emission probabilities')
        assert_almost_equal(np.array(m_pi), np.array(pi), decimal=decimal, err_msg='Wrong initial probabilities')
