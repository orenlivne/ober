'''
============================================================
Test the numpy hashable wrapper.

Created on September 14, 2012
@author: Oren Livne <livne@uchicago.edu>
@see http://machineawakening.blogspot.com/2011/03/making-numpy-ndarrays-hashable.html
============================================================
'''
import unittest, numpy as np
from numpy.ma.testutils import assert_equal
from hashable import Hashable
from numpy import random

class TestHashable(unittest.TestCase):
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_dict_of_hashable(self):
        '''Test creating a dictionary with a hashable entry.'''
        b = Hashable(np.arange(0, 2024))
        d = {}
        d[b] = 'bar'
        assert_equal(d[b], 'bar', 'Wrong dictionary entry')

    def test_hashes_of_equal_objects_are_equal(self):
        '''Test reproducability of hashing long arrays.'''
        a = Hashable(np.arange(0, 10000))
        b = Hashable(np.arange(0, 10000))
        assert_equal(hash(a), hash(b), 'Hash should be reproducible')

    def test_hashes_of_equal_objects_from_slice_are_equal(self):
        '''Test reproducability of hashing boolean rows of a matrix obtained from slicing.'''
        #h = np.mod(np.arange(0, 100).reshape(100, 3), 2)
        h = random.randint(0, 2, (4, 100, 2))
        h = h[:, :, 1].transpose()
        h = h[:, np.arange(0, 3)]
        h = (h == 1)
        def get_index(h):
            for j in xrange(1, h.shape[0]):
                if np.all(h[j] == h[0]):
                    return j
            return None
        j = get_index(h)
        assert_equal(h[0], h[j], 'Test objects should be equal, if not, rerun test')
        #def hash_array(x): return sha1(buffer(x.view(np.uint8))).hexdigest()
        def hash_array(x): x.tostring()
        assert_equal(hash_array(h[0]), hash_array(h[j]), 'Hash should be reproducible')

    def test_hash_rows(self):
        '''Test hashing the rows of a numpy matrix.'''
        a = np.array([[1, 2, 3], [4, 5, 6], [1, 2, 3], [7, 8, 9]], dtype=int)
        d = dict(((Hashable(x), i) for (i, x) in enumerate(a)))
        assert_equal(len(d), 3, 'Wrong hash table size')
        # Hash value should be the last-encountered index of the corresponding row 
        assert_equal(d[Hashable(np.array([1, 2, 3]))], 2, 'Wrong hash table entry')
