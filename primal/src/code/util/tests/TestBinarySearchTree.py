'''
============================================================
Test binary search trees. 

Created on July 17, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''

import unittest
import numpy as np
from numpy.ma.testutils import assert_equal
from bst import BinarySearchTree
from util import alternating_order, optimal_insertion_order, nearest_neighbor,\
    sequence_to_tree, nearest_neighbor_multiple
from numpy.ma.core import sort

class TestBinarySearchTree(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
 
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_queries(self):
        '''Create a simple test tree and test BST queries.''' 
        b = TestBinarySearchTree._simple_tree()
        assert_equal(b.minimum(), 1, 'Wrong BST min value search')
        assert_equal(b.maximum(), 14, 'Wrong BST min value search')
        
        assert_equal(b.find(3), 3, 'Wrong BST search')
        # Find inexistent element
        try:
            b.find(-1)
            raise ValueError('Wrong BST search')
        except KeyError:
            # OK
            pass
        
    def test_find_largest_le(self):
        '''Test our additional search function.'''
        b = TestBinarySearchTree._simple_tree()
        assert_equal(b.find_largest_le(1.5), 1, 'Wrong BST largest element less than search')
        assert_equal(b.find_largest_le(6.5), 6, 'Wrong BST largest element less than search')
        assert_equal(b.find_largest_le(8.5), 8, 'Wrong BST largest element less than search')
        assert_equal(b.find_largest_le(14.5), 14, 'Wrong BST largest element less than search')
        
        assert_equal(b.find_largest_le(1), 1, 'Wrong BST largest element less than search')
        assert_equal(b.find_largest_le(7), 7, 'Wrong BST largest element less than search')
        assert_equal(b.find_largest_le(14), 14, 'Wrong BST largest element less than search')

        assert_equal(b.find_largest_le(0), None, 'Wrong BST largest element less than search')
        
    def test_find_smallest_ge(self):
        '''Test our additional search function.'''
        b = TestBinarySearchTree._simple_tree()
        assert_equal(b.find_smallest_ge(0), 1, 'Wrong BST smallest element greater than search')
        assert_equal(b.find_smallest_ge(1.5), 3, 'Wrong BST smallest element greater than search')
        assert_equal(b.find_smallest_ge(6.5), 7, 'Wrong BST smallest element greater than search')
        assert_equal(b.find_smallest_ge(8.5), 10, 'Wrong BST smallest element greater than search')
        
        assert_equal(b.find_smallest_ge(1), 1, 'Wrong BST smallest element greater than search')
        assert_equal(b.find_smallest_ge(7), 7, 'Wrong BST smallest element greater than search')
        assert_equal(b.find_smallest_ge(14), 14, 'Wrong BST smallest element greater than search')

        assert_equal(b.find_smallest_ge(14.5), None, 'Wrong BST smallest element greater than search')
    
    def test_alternating_order(self):
        '''Test alternating order generation.'''
        assert_equal(alternating_order(5), [2, 1, 3, 0, 4], 'Wrong alternating order for odd n')
        assert_equal(alternating_order(6), [3, 2, 4, 1, 5, 0], 'Wrong alternating order for even n')
        assert_equal(alternating_order(1), [0], 'Wrong alternating order for n=1')
        assert_equal(alternating_order(0), [], 'Wrong alternating order for n=0')

    def test_alternating_order_insertion(self):
        '''Test alternating order generation.'''
        
        # Minimum possible depth = [log2(n)]+1 where n=size of data, so here 4
        data = np.array([5, 1, 2, 7, 8, 3, 4, 6, 9])
        n = np.size(data)
        
        # Original list order
        self._test_tree_depth(data, 5)
        # Ordered insertion is worst case!
        sorted_data = sort(data)
        self._test_tree_depth(sorted_data, 9)
        # Alternating order is not so good either 
        self._test_tree_depth(sorted_data[alternating_order(n)], 5)
        # Neither does random order guarantee anything
        self._test_tree_depth(data[np.random.permutation(n)], 4, 9)

        # Best order. O(n log n) complexity to sort the list.
        self._test_tree_depth(sorted_data[optimal_insertion_order(n)], 4)
        self._test_tree_depth([8, 4, 2, 1, 3, 6, 5, 7, 9], 4)
        
        # Best order holds for a general list. All numbers in the list must be different.
        n = 100
        data = np.random.uniform(size=n)#np.random.permutation(n)
        sorted_data = sort(data)
        depth = int(np.floor(np.log2(n)))+1
        self._test_tree_depth(sorted_data[optimal_insertion_order(n)], depth)
        # Same code using a util function
        b = sequence_to_tree(sorted_data)
        assert_equal(b.depth(), depth, 'Wrong tree depth')
    
    def test_nearest_neighbor(self):
        '''Test our additional search function.'''
        locations = np.array([0, 5, 6.5, 8, 11])
        assert_equal(0, nearest_neighbor(-0.5, locations), 'Wrong nearest neighbor index')
        assert_equal(0, nearest_neighbor(2.4, locations), 'Wrong nearest neighbor index')
        assert_equal(1, nearest_neighbor(3, locations), 'Wrong nearest neighbor index')
        assert_equal(3, nearest_neighbor(7.5, locations), 'Wrong nearest neighbor index')
        assert_equal(1, nearest_neighbor(5.0, locations), 'Wrong nearest neighbor index')
        assert_equal(4, nearest_neighbor(10, locations), 'Wrong nearest neighbor index')
        assert_equal(4, nearest_neighbor(11.5, locations), 'Wrong nearest neighbor index')
        assert_equal([0, 3], list(nearest_neighbor_multiple([0, 7.5], locations)), 'Wrong nearest neighbor index')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    @staticmethod
    def _simple_tree():
        '''Create a simple test tree. See http://en.wikipedia.org/wiki/Binary_search_tree'''
        return BinarySearchTree(values=[8, 3, 1, 6, 4, 7, 10, 14, 13])
        
    def _test_tree_depth(self, data, expected_depth, upper_bound=None):
        '''Test tree depth for a certain insertion order (data).'''
        b = BinarySearchTree(values=data)
        #print data
        #print b.pprint()
        if upper_bound is None:
            assert_equal(b.depth(), expected_depth, 'Wrong tree depth')
        else:
            d = b.depth()
            self.assertTrue(d >= expected_depth and d <= upper_bound, 'Tree depth %d out of expected range [%d,%d]' % (d, expected_depth, upper_bound,))
