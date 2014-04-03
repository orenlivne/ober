'''
============================================================
Test the GHMM Hidden Markov Model (HMM) library.

Created on November 20, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, statutil, hmm
from numpy.testing.utils import assert_almost_equal
from numpy.ma.testutils import assert_equal

class TestHmm(unittest.TestCase):
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
        A = np.array([ [.5, .5], [.5, .5]])
        # Emission probabilities
        B = np.array([ [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, ], \
                         [ 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 / 2 ] ])
        # Symbols
        V = [1, 2, 3, 4, 5, 6]

        # Model
        m = hmm.HMM(2, A=A, B=B, V=V)
        TestHmm.assert_model_matrices_almost_equal(m, (A, B, [0.5, 0.5])) 
                                                   
    def test_dishonest_casino(self):
        '''Dishonest Casino Example.'''
        # Create transition probability matrix
        A = np.array([[0.99, 0.01],
                      [0.01, 0.99]])
        # Create observable probability distribution matrix. Casino biased toward "6" in state "1".        
        B = statutil.scale_row_sums(np.array([[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ],
                                              [ 1.0, 1.0, 1.0, 1.0, 1.0, 5.0 ]]))
        # Create set of all observable symbols
        V = [1, 2, 3, 4, 5, 6]
    
        # Instantiate an HMM, note Pi is uniform probability distribution by default
        m = hmm.HMM(2, A=A, B=B, V=V)
        
        Obs = [ 1, 2, 3, 4, 5, 2, 1, 6, 6, 6, 5, 6 ]
        log_prob_Obs, Alpha, c = hmm.forward(m, Obs, scaling=1)
        assert_almost_equal(log_prob_Obs, -20.9468006, decimal=5, err_msg='Wrong observation probability')
        
        Q_star, _, _ = hmm.viterbi(m, Obs, scaling=1)
        assert_equal(Q_star, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 'Wrong Viterbi path')

        Beta = hmm.backward(m, Obs, c)
        Gamma, Q_star = hmm.individually_optimal_states(Alpha, Beta)
        assert_almost_equal(Gamma,
                            [[0.63711364302936, 0.6348934929050587, 0.6271179131667495, 0.6117100305977996, 0.5845543683193845, 0.5383975935172204, 0.46091113744414974, 0.3313982095474306, 0.28864618346708165, 0.27562909135388625, 0.27498372625848855, 0.26932891011973825], [0.36288635697064003, 0.3651065070949412, 0.3728820868332506, 0.38828996940220045, 0.4154456316806155, 0.4616024064827796, 0.5390888625558502, 0.6686017904525694, 0.7113538165329184, 0.7243709086461138, 0.7250162737415115, 0.7306710898802617]],
                            decimal=5, err_msg='Wrong state probabilities')        
        assert_equal(Q_star, [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1], 'Wrong individually-optimal states')
        
    def test_dishonest_casino_larger_transition_p(self):
        '''Dishonest Casino Example.'''
        # Create transition probability matrix
        A = np.array([[0.9, 0.1],
                      [0.1, 0.9]])
        # Create observable probability distribution matrix. Casino biased toward "6" in state "1"
        B = statutil.scale_row_sums(np.array([[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ],
                                              [ 1.0, 1.0, 1.0, 1.0, 1.0, 5.0 ]]))
        # Create set of all observable symbols
        V = [1, 2, 3, 4, 5, 6]
    
        # Instantiate an HMM, note Pi is uniform probability distribution by default
        m = hmm.HMM(2, A=A, B=B, V=V)
        
        Obs = [ 1, 2, 3, 4, 5, 2, 1, 6, 6, 6, 5, 6 ]
        log_prob_Obs, Alpha, c = hmm.forward(m, Obs, scaling=1)
        assert_almost_equal(log_prob_Obs, -20.124, decimal=3, err_msg='Wrong observation probability')
        
        Q_star, _, _ = hmm.viterbi(m, Obs, scaling=1)
        assert_equal(Q_star, [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1], err_msg='Wrong Viterbi path')

        Beta = hmm.backward(m, Obs, c)
        Gamma, Q_star = hmm.individually_optimal_states(Alpha, Beta)
        assert_almost_equal(Gamma,
                            [[0.8189770516168013, 0.8482906260695058, 0.8525027084764197, 0.8329611652077556, 0.7834127024175411, 0.6880018120129073, 0.5161970090643716, 0.2130207566284025, 0.12024202874950358, 0.10797060639721641, 0.15902649827833876, 0.14930464162738483], [0.18102294838319855, 0.15170937393049422, 0.14749729152358024, 0.16703883479224435, 0.21658729758245884, 0.31199818798709256, 0.4838029909356284, 0.7869792433715975, 0.8797579712504964, 0.8920293936027837, 0.8409735017216613, 0.8506953583726152]],
                            decimal=5, err_msg='Wrong state probabilities')        
        assert_equal(Q_star, [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1], 'Wrong individually-optimal states')
        
    def test_train_model(self):
        '''Dishonest Casino Example - EM algorithm.'''
        # Create transition probability matrix
        A = np.array([[0.99, 0.01],
                      [0.01, 0.99]])
        # Create observable probability distribution matrix. Casino biased toward "6" in state "1".        
        B = statutil.scale_row_sums(np.array([[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ],
                                              [ 1.0, 1.0, 1.0, 1.0, 1.0, 5.0 ]]))
        # Create set of all observable symbols
        V = [1, 2, 3, 4, 5, 6]
    
        # Instantiate an HMM, note Pi is uniform probability distribution by default
        m = hmm.HMM(2, A=A, B=B, V=V)
        
        Obs = [ 1, 2, 3, 4, 5, 2, 1, 6, 6, 6, 5, 6 ]
        c = [Obs]
        hmm.baum_welch(m, c, epochs=15, graph=False)
        TestHmm.assert_model_matrices_almost_equal(m, 
                                                   ([[0.856658708052639, 0.14334129194736125], [2.454940916925095e-16, 1.0]],
                                                    [[0.28329354031233306, 0.2866825838637413, 0.14334129194736112, 0.14334129194736112, 0.14334129192821368, 9.896623857864685e-13], [0.004706380704415612, 4.3023359620169447e-11, 3.2510873580469717e-111, 1.2201233032249015e-54, 0.19905872387205914, 0.7962348953805019]],
                                                    [1.0, 4.364785210913299e-122]))
                
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    @staticmethod
    def assert_model_matrices_almost_equal(m, m_expected, decimal=3):
        '''Check that a model''s matrices equal a set of expected matrices.'''
        (A, B, Pi) = m_expected
        assert_almost_equal(m.A, np.array(A), decimal=decimal, err_msg='Wrong transition probabilities')
        assert_almost_equal(m.B, np.array(B), decimal=decimal, err_msg='Wrong emission probabilities')
        assert_almost_equal(m.Pi, np.array(Pi), decimal=decimal, err_msg='Wrong initial probabilities')
