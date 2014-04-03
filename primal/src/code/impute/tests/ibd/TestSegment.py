'''
============================================================
Test basic discrete interval operations. 

Created on August 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, impute as im
from numpy.ma.testutils import assert_equal

class TestSegment(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------

    def setUp(self):
        '''Load single nuclear family test case.'''
        pass
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_group_to_disjoint_scalar(self):
        '''Check that breaking intervals of pairwise equality between values into disjoint intervals
        with globally-equal-value sets works. Values = scalars'''
        segments = [((0, 50), (1, 2)),
                    ((50, 100), (1, 3)),
                    ((50, 300), (2, 4)),
                    ((100, 200), (3, 4))]
        broken = [((0, 50), set([1, 2])),
                    ((50, 100), set([1, 3])),
                    ((50, 100), set([2, 4])),
                    ((100, 200), set([2, 3, 4])),
                    ((200, 300), set([2, 4]))]
        assert_equal(list(im.segment.group_to_disjoint(segments)), broken, 'Wrong segment breaking')
        
    def test_group_to_disjoint_tuple(self):
        '''Check that breaking intervals of pairwise equality between values into disjoint intervals
        with globally-equal-value sets works. Values = tuples.'''
        segments = [((0, 50), ((1, 0), (2, 0))),
                    ((50, 100), ((1, 0), (3, 1))),
                    ((50, 300), ((2, 0), (4, 1))),
                    ((100, 200), ((3, 1), (4, 1)))]
        broken = [((0, 50), frozenset([(2, 0), (1, 0)])),
                    ((50, 100), frozenset([(2, 0), (4, 1)])),
                    ((50, 100), frozenset([(1, 0), (3, 1)])),
                    ((100, 200), frozenset([(2, 0), (3, 1), (4, 1)])),
                    ((200, 300), frozenset([(2, 0), (4, 1)]))]
        assert_equal(list(im.segment.group_to_disjoint(segments)), broken, 'Wrong segment breaking')

    def test_group_to_disjoint_merge(self):
        '''Check that we correctly merge identical consecutive segments.'''
        segments = [((0, 50), (1, 2)),
                    ((50, 300), (2, 3)),
                    ((50, 100), (3, 4)),
                    ((100, 200), (2, 4)),
                    ((200, 300), (3, 4)),
                    ((300, 400), (1, 4))]
        broken = [((0, 50), set([1, 2])),
                    ((50, 300), set([2, 3, 4])),
                    ((300, 400), set([1, 4]))]
        assert_equal(list(im.segment.group_to_disjoint(segments, merge_consecutive=True)), broken, 'Wrong segment breaking')
   
   
    def test_group_to_disjoint_no_merge(self):
        '''Check that we correctly merge identical consecutive segments.'''
        segments = [((0, 50), (1, 2)),
                    ((50, 300), (2, 3)),
                    ((50, 100), (3, 4)),
                    ((100, 200), (2, 4)),
                    ((200, 300), (3, 4)),
                    ((300, 400), (1, 4))]
        broken = [((0, 50), set([1, 2])),
                    ((50, 100), set([2, 3, 4])),
                    ((100, 200), set([2, 3, 4])),
                    ((200, 300), set([2, 3, 4])),
                    ((300, 400), set([1, 4]))]
        assert_equal(list(im.segment.group_to_disjoint(segments)), broken, 'Wrong segment breaking')

    def test_segments_with_value(self):
        '''Check that we correctly identify IBS segments (0-segments) within a binary array.'''
        A = np.array([1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1], dtype=np.byte)
        assert_equal(im.segment.segments_with_value(A, 0, 0), [(1, 3), (4, 7), (9, 13)], 'Wrong zero segments')
        assert_equal(im.segment.segments_with_value(A, 0, 3), [(4, 7), (9, 13)], 'Wrong zero segments - minimum length size')
        assert_equal(im.segment.segments_with_value(A, 0, 5), [], 'Wrong zero segments - no segments found')

        A = np.array([1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0], dtype=np.byte)
        assert_equal(im.segment.segments_with_value(A, 0, 0), [(1, 3), (4, 7), (9, 16)], 'Wrong zero segments - right boundary case')

        A = np.array([0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1], dtype=np.byte)
        assert_equal(im.segment.segments_with_value(A, 0, 0), [(0, 3), (4, 7), (9, 13)], 'Wrong zero segments - left boundary case')
        assert_equal(im.segment.segments_with_value(A, 1, 3), [(13, 16)], 'Wrong zero segments - value != 0')

        A = np.array(np.zeros(10), dtype=np.byte)
        assert_equal(im.segment.segments_with_value(A, 0, 0), [(0, 10)], 'Wrong zero segments - zero array')

    def test_segment_set_op(self):
        '''Test segment set operations.''' 
        a = [(0, 2), (4, 7), (10, 100)]
        b = [(1, 5), (6, 20)]

        assert_equal(im.segment.segment_set_op(a, a, 'and'), a, 'Wrong segment set intersection')
        assert_equal(im.segment.segment_set_op(a, a, 'or'), a, 'Wrong segment set intersection')

        assert_equal(im.segment.segment_set_op([], [], 'or'), [], 'Wrong segment set intersection')
        assert_equal(im.segment.segment_set_op([], [], 'and'), [], 'Wrong segment set intersection')
        assert_equal(im.segment.segment_set_op(a, [], 'or'), a, 'Wrong segment set intersection')
        assert_equal(im.segment.segment_set_op(a, [], 'and'), [], 'Wrong segment set intersection')
        assert_equal(im.segment.segment_set_op([], a, 'or'), a, 'Wrong segment set intersection')
        assert_equal(im.segment.segment_set_op([], a, 'and'), [], 'Wrong segment set intersection')

        assert_equal(im.segment.segment_set_op(a, b, 'and'), [(1, 2), (4, 5), (6, 7), (10, 20)], 'Wrong segment set intersection')
        assert_equal(im.segment.segment_set_op(a, b, 'and'), [(1, 2), (4, 5), (6, 7), (10, 20)], 'Wrong segment set intersection')
        assert_equal(im.segment.segment_set_op(b, a, 'or'), [(0, 100)], 'Wrong segment set union')
        assert_equal(im.segment.segment_set_op(b, a, 'or'), [(0, 100)], 'Wrong segment set union')
        
    def test_segment_union(self):
        '''Test unionning segments.''' 
        assert_equal(im.segment.union([(0, 2), (4, 7), (10, 100)]), [(0, 2), (4, 7), (10, 100)], 'Wrong segment union')
        assert_equal(im.segment.union([(0, 7), (4, 7), (10, 100)]), [(0, 7), (10, 100)], 'Wrong segment union')
        assert_equal(im.segment.union([(0, 7), (4, 6), (2, 5), (10, 100)]), [(0, 7), (10, 100)], 'Wrong segment union')
        assert_equal(im.segment.union([(0, 7), (4, 9), (2, 5), (10, 100)]), [(0, 9), (10, 100)], 'Wrong segment union')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
