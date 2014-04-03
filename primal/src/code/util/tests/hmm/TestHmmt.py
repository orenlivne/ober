'''
============================================================
Test the time-dependent Hidden Markov Model (HMM) library.

Created on November 20, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, hmmt, statutil
from numpy.testing.utils import assert_almost_equal
from numpy.ma.testutils import assert_equal

class TestHmmt(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Create a simple dice rolling HMM'''
        # adjust the precision of printing float values
        np.set_printoptions(precision=4)
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_create_model(self):
        '''Based on Mike's DC example.'''
        # Transition probabilities 
        c = 1.0
        A = lambda t: np.array([[.5 + c / (t + 10), .5 - c / (t + 10)],
                                [.5 - c / (t + 10), .5 + c / (t + 10)]])
        # Emission probabilities (actually need to sum up to 1 to be proper probabilities! Doesn't matter here)
        B = lambda t: np.array([[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, ],
                                [ 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 / 2 ]])
        # Symbols
        V = [1, 2, 3, 4, 5, 6]

        # Model
        m = hmmt.HMM(2, A=A, B=B, V=V)
        TestHmmt.assert_model_matrices_almost_equal(m, 0, ([[0.6, 0.4], [0.4, 0.6]],
                                                        [[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, ],
                                                         [ 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 / 2 ]],
                                                        [0.5, 0.5]))
        
    def test_dishonest_casino(self):
        '''Dishonest Casino Example.'''
        # Create set of all observable symbols
        V = [1, 2, 3, 4, 5, 6]
    
        # Instantiate an HMM, note Pi is uniform probability distribution by default
        m = hmmt.HMM(2, A=A, B=B, V=V)
        
        Obs = [ 1, 2, 3, 4, 5, 2, 1, 6, 6, 6, 5, 6 ]
        log_prob_Obs, Alpha, c = hmmt.forward(m, Obs, scaling=1)
        assert_almost_equal(log_prob_Obs, -20.918331960920177, decimal=5, err_msg='Wrong observation probability')
        
        Q_star, _, _ = hmmt.viterbi(m, Obs, scaling=1)
        assert_equal(Q_star, np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]), 'Wrong Viterbi path')

        Beta = hmmt.backward(m, Obs, c)
        Gamma, Q_star = hmmt.individually_optimal_states(Alpha, Beta)
        assert_almost_equal(Gamma,
                            [[0.6497764611558418, 0.6476464636785156, 0.6395341323290045, 0.6232731859090421, 0.594520672835795, 0.5455970624010458, 0.4634351125974832, 0.3259525777265338, 0.28060687540229046, 0.2668719322118299, 0.26641411220848843, 0.26052432611296034], [0.35022353884415813, 0.3523535363214843, 0.3604658676709956, 0.3767268140909578, 0.405479327164205, 0.4544029375989542, 0.5365648874025167, 0.6740474222734663, 0.7193931245977097, 0.7331280677881701, 0.7335858877915117, 0.7394756738870396]],
                            decimal=5, err_msg='Wrong state probabilities')        
        assert_equal(Q_star, [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1], 'Wrong individually-optimal states')
        
    def test_dishonest_casino_larger_transition_p(self):
        '''Dishonest Casino Example.'''
        # Create transition probability matrix
        # Create set of all observable symbols
        V = [1, 2, 3, 4, 5, 6]
    
        # Instantiate an HMM, note Pi is uniform probability distribution by default
        m = hmmt.HMM(2, A=A2, B=B, V=V)
        
        Obs = [ 1, 2, 3, 4, 5, 2, 1, 6, 6, 6, 5, 6 ]
        log_prob_Obs, Alpha, c = hmmt.forward(m, Obs, scaling=1)
        assert_almost_equal(log_prob_Obs, -20.12266, decimal=3, err_msg='Wrong observation probability')
        
        Q_star, _, _ = hmmt.viterbi(m, Obs, scaling=1)
        assert_equal(Q_star, [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1], err_msg='Wrong Viterbi path')

        Beta = hmmt.backward(m, Obs, c)
        Gamma, Q_star = hmmt.individually_optimal_states(Alpha, Beta)
        assert_almost_equal(Gamma,
                            [[0.8187172809335708, 0.8482909847967097, 0.8526704789882301, 0.8332628435868141, 0.7838346435680826, 0.6885099864086844, 0.5166714952004011, 0.21280359871744636, 0.11993462067436834, 0.10785580156479256, 0.15963081036800023, 0.1497144396249583], 
                             [0.1812827190664292, 0.1517090152032903, 0.14732952101176994, 0.16673715641318582, 0.21616535643191742, 0.31149001359131556, 0.4833285047995989, 0.7871964012825536, 0.8800653793256317, 0.8921441984352074, 0.8403691896319997, 0.8502855603750418]],
                            decimal=5, err_msg='Wrong state probabilities')
        assert_almost_equal(np.sum(Gamma[np.array([0,1]), :], axis=0), np.ones(len(Obs)),
                                decimal=3, err_msg='Gammas must sum up to 1 at each time')      
        assert_equal(Q_star, [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1], 'Wrong individually-optimal states')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    @staticmethod
    def assert_model_matrices_almost_equal(m, t, m_expected, decimal=3):
        '''Check that a model''s matrices equal a set of expected matrices.'''
        (A, B, Pi) = m_expected
        assert_almost_equal(m.A(t), np.array(A), decimal=decimal, err_msg='Wrong transition probabilities')
        assert_almost_equal(m.B(t), np.array(B), decimal=decimal, err_msg='Wrong emission probabilities')
        assert_almost_equal(m.Pi, np.array(Pi), decimal=decimal, err_msg='Wrong initial probabilities')

D = 1000.0

def A(t): 
    return np.array([[0.99 - 1.0 / (t + D), 0.01 + 1.0 / (t + D)],
                     [0.01 + 1.0 / (t + D), 0.99 - 1.0 / (t + D)]])

def A2(t): 
    return np.array([[0.9 - 1.0 / (t + D), 0.1 + 1.0 / (t + D)],
                     [0.1 + 1.0 / (t + D), 0.9 - 1.0 / (t + D)]])

def B(t):    
    # 1.0reate observable probability distribution matrix. 1.0asino biased toward "6" in state "1".        
    return statutil.scale_row_sums(np.array([[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ],
                                             [ 1.0 - 1.0 / (t + D), 1.0, 1.0, 1.0, 1.0, 5.0 + 1.0 / (t + D) ]]))
