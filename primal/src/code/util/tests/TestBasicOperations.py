'''
============================================================
Test basic python behaviour.

Created on May 30, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, util, unittest
from tests.SimpleClasses import Clazz1
from numpy.testing.utils import assert_equal
from util import dict_invert, max_except_center_filter
from scipy.ndimage.filters import generic_filter

class TestBasicOperations(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------   
    
    def test_arg_pass_by_reference(self):
        '''Test whether an method object argument is passed by reference.'''
        c = Clazz1()
        assert c.x == 200, 'Property not set properly'
        TestBasicOperations.__my_method(c)
        assert c.x == 300, 'Property not updated after returning from method'
        
    def test_numpy_vectorize(self):
        '''Test elementwise-applying a dictionary to a numpy array.'''
        d = dict(zip(range(0, 4), range(4, 8)))
        a = np.array([[0, 1], [2, 3]])
        f = np.vectorize(d.__getitem__)
        assert_equal(f(a), [[4, 5], [6, 7]], 'Wrong vectorize result')

    def test_dict_invert(self):
        '''Test inverting a 1:1 dictionary.'''
        d = dict(zip(range(0, 4), range(4, 8)))
        d_inv = dict(zip(range(4, 8), range(0, 4)))
        assert_equal(dict_invert(d), d_inv, 'Wrong dictionary inversion')

    def test_generic_filter(self):
        '''Custom digitla filter.'''
        filtered_array = generic_filter(np.array([0., 1, 2, 3, 4, 5, 6, 7]), my_function, 3)
        assert_equal(filtered_array, [0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.25], 'Wrong custom filter result')

    def test_filter_max_except_center(self):
        '''Custom digitla filter.'''
        assert_equal(max_except_center_filter(np.array([0, 0, 1, 0, 0, 1]), 3), [0, 0, 0, 0, 0, 0], 'Wrong custom filter result')
        assert_equal(max_except_center_filter(np.array([0, 0, 1, 1, 0, 1]), 5), [0, 0, 1, 1, 0, 1], 'Wrong custom filter result')
    
    def test_hamming_distance(self):
        '''Test the bitwise implementation of the hamming distance function.'''
        num_experiments = 100
        a = np.reshape(np.random.random_integers(0, 255, size=2 * num_experiments),
                       (num_experiments, 2)).astype(np.uint8)
        for row in a:
            (i, j) = row[0:2]
            assert_equal(util.hamming_distance(i, j), TestBasicOperations.__hamming_lazy(i, j, 8),
                         'Wrong hamming distance computation for i=%d, j=%d' % (i, j))
        
    def test_nested_list_comprehension(self):
        '''Test nested list comprehension expressions.'''
        d = {1: [2, 3], 2: [4, 5]}
        assert_equal([(k, x) for (k, v) in d.iteritems() for x in v], [(1, 2), (1, 3), (2, 4), (2, 5)],
                     'Wrong list comprehension result')
         
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    @staticmethod
    def __my_method(c):
        c.x = 300
        return c

    @staticmethod
    def __hamming_lazy(a, b, bits=32): 
        x = (a ^ b) & ((1 << bits) - 1) 
        tot = 0 
        while x:
            tot += x & 1 
            x >>= 1 
        return tot

def my_function(window):
    return (window[0] + window[2]) / 4.0

