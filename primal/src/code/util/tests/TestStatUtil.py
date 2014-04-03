'''
============================================================
Test our statistical utilities.

Created on October 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, statutil
from numpy.testing.utils import assert_equal
from scipy.stats.stats import linregress
from numpy.ma.testutils import assert_almost_equal

class TestStatUtil(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    SEED = 0
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Generate test data.'''
        np.random.seed(TestStatUtil.SEED)
            
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------       
    def test_multinomial_elementwise_vector(self):
        '''Test creating multinomial variables (r=1).'''
        (m, n) = (20, 5)
        p = statutil.random_row_stochastic((m, n))
        x = statutil.multinomial_elementwise(p)
        assert_equal(x, [3, 1, 3, 4, 0, 2, 2, 3, 1, 4, 2, 3, 2, 1, 4, 2, 4, 3, 4, 3],
                     'Wrong random multinomial generated')

    def test_multinomial_elementwise_matrix(self):
        '''Test creating multinomial variables (r > 1).'''
        (m, n) = (20, 5)
        p = statutil.random_row_stochastic((m, n))
        x = statutil.multinomial_elementwise(p, 2)
        assert_equal(x, [[3, 3], [1, 2], [3, 3], [4, 3], [0, 1], [2, 2], [2, 0], [3, 1], [1, 2],
                         [4, 0], [2, 2], [3, 1], [2, 1], [1, 1], [4, 2], [2, 3], [4, 3], [3, 4],
                         [4, 3], [3, 1]],
                     'Wrong random multinomial generated')
    
    def test_multinomial_elementwise_distribution(self):
        '''Verify that the created variables approach a multinomial distribution for large numbers
        of samples.'''
        (m, n, k) = (6, 5, 1)
        r = 2 ** np.arange(4, 17)
        p = statutil.random_row_stochastic((m, n))
        #p = statutil.scale_row_sums(np.ones((m, n)))
        error = np.zeros((len(r),))
        for (i, r_val) in enumerate(r):
            for _ in xrange(k):
                x = statutil.multinomial_elementwise(p, r_val)
                # Root-mean-square-error of observed frequencies w.r.t. desired frequencies
                error[i] += statutil.norm_frobenius_scaled(statutil.hist(x, n) / (1.0 * r_val) - p)
            error[i] /= (1.0 * k)
        # Validate the model error of the central limit theorem: C*r^(-0.5).
        # This is a consequence of the Central Limit Theorem. We are making k experiments for
        # each value of n. Even if k=1, there's a 95% chance that we are within ~1.6 standard deviations
        # from the mean of the normal distribution sqrt(n)*[observed freq variable - p[i,j]] for each
        # entry j of a row i of the matrix p. So if row i's stddev is s[i], the sum of squared errors
        # should be (with 95% confidence) <= n * (1.96*s[i])^2. So 
        # C <= sqrt(sum(n * (1.5*s[i])^2)_i / (m*n)) = 1.96 * sqrt(s[i]^2/m).
        # See http://en.wikipedia.org/wiki/Central_limit_theorem
        alpha, c, r_value, _, _ = linregress(np.log(r), np.log(error))
        c = np.exp(c)
#        print c , 1.96 * np.linalg.linalg.norm(np.sum(p * np.arange(p.shape[1]) ** 2, axis=1) - 
#                                                          np.sum(p * np.arange(p.shape[1]), axis=1) ** 2,
#                                                          2) / np.sqrt(p.shape[0]),
        assert_almost_equal(alpha, -0.5, decimal=1, err_msg='Unexpected error term growth power')
        self.assertTrue(c <= 1.96 * np.linalg.linalg.norm(np.sum(p * np.arange(p.shape[1]) ** 2, axis=1) - 
                                                          np.sum(p * np.arange(p.shape[1]), axis=1) ** 2,
                                                          2) / np.sqrt(p.shape[0]),
                        'Error term coefficient outside 95% confidence interval')
        self.assertTrue(abs(r_value) > 0.99, 'Error does not fit a power law in sample size')
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __test_segment_range(self, n, step, c):
        '''Test getting the start and end of equidistant intervals in an item collection.'''
        pass
