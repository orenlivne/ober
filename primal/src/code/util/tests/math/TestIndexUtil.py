'''
============================================================
Test array index Cython utilities.

Created on December 4, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, timeit, os
from numpy.ma.testutils import assert_equal
from utilities.math.index_util import first_index_long, first_occurrence_index_long, first_occurrence_index_byte
from utilities.math import index_util
from utilities.math.filter import median_filter

class TestIndexUtil(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------

    def setUp(self):
        '''Set up paths.'''
        self.dir = os.path.dirname(__file__)
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_first_index(self):
        '''Test finding the index of the first occurrence of a value in an array.'''
        a = np.random.random_integers(0, 10000, 1000000)
        value = 57
        index = first_index_long(a, value)
        assert_equal(a[index], value, 'Wrong index finding')

    def test_time_first_index(self):
        t1 = timeit.Timer('np.where(a==57)[0][0]', 'import numpy as np; np.random.seed(1); a = np.random.random_integers(0, 100000, 1000000)')
        t1 = t1.timeit(100) / 100
        t2 = timeit.Timer('lib.find_index(57, a.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), len(a))', 'import numpy as np; np.random.seed(1); a = np.random.random_integers(0, 100000, 1000000); import ctypes; lib = ctypes.CDLL("%s/find_index.o.dylib"); lib.find_index.restype = ctypes.c_long; lib.find_index.argtypes = (ctypes.c_long, ctypes.POINTER(ctypes.c_long), ctypes.c_long) ' % index_util.DIR)
        t2 = t2.timeit(100) / 100
        self.assertTrue(t2 < 0.1 * t1, 'Speed up of at least x10 is expected with Cython find_index() call')

    def test_first_occurrence_index_long(self):
        '''Test finding the index of the first occurrence of a value in an array.'''
        a = np.array([0, 1, 2, 3, 4, 5, 6, 3, 2, 1])
        assert_equal(first_occurrence_index_long(a, 3, 0, 1), 3, 'Wrong index finding, forward')
        assert_equal(first_occurrence_index_long(a, 3, len(a) - 1, -1), 7, 'Wrong index finding, backward')
        assert_equal(first_occurrence_index_long(a, 3, 5, -1), 3, 'Wrong index finding, from middle backward')
        assert_equal(first_occurrence_index_long(a, 3, 5, 1), 7, 'Wrong index finding, from middle forward')
        assert_equal(first_occurrence_index_long(a, 10, 0, 1), -1, 'Wrong index finding, from middle forward, value not in array')

    def test_first_occurrence_index_int8(self):
        '''Test finding the index of the first occurrence of a value in an array.'''
        a = np.array([0, 0, 0, 1, 0, 0, 0, 1, 0, 0]).astype('int8')
        assert_equal(first_occurrence_index_byte(a, 1, 0, 1), 3, 'Wrong index finding, forward')
        assert_equal(first_occurrence_index_byte(a, 1, len(a) - 1, -1), 7, 'Wrong index finding, backward')
        assert_equal(first_occurrence_index_byte(a, 1, 5, -1), 3, 'Wrong index finding, from middle backward')
        assert_equal(first_occurrence_index_byte(a, 1, 5, 1), 7, 'Wrong index finding, from middle forward')
        assert_equal(first_occurrence_index_byte(a, 10, 0, 1), -1, 'Wrong index finding, from middle forward, value not in array')

    def test_median_filter(self):
        '''Test finding the index of the first occurrence of a value in an array.'''
        a = np.ones((2, 3), dtype=np.uint8)
        print a
        median_filter(a, 3)
        print a
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
