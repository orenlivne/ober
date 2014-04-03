'''
============================================================
Test array neighbor index finding routines. 

Created on May 30, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, util, unittest
from numpy.ma.testutils import assert_equal

class TestUtilNeighbor(unittest.TestCase):
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_nearest_neighbor_in_sorted_arrays(self):
        a = np.sort(np.unique(np.random.randint(1000, size=120)))
        b = np.sort(np.unique(np.random.randint(1000, size=150)))
        actual = list(util.nearest_neighbor_in_sorted_arrays(a, b))
        expected = [np.abs(a - x).argmin() for x in b]
        assert_equal(actual, expected, 'Wrong nearest neighbor location in sorted arrays')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
