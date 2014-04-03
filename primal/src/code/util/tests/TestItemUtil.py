'''
============================================================
Test basic python behaviour.

Created on May 30, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, test_util as tu, itemutil as iu, numpy as np
from numpy.testing.utils import assert_equal

class TestItemUtil(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    __FILE = tu.abs_path('misc/hutt.chr22.snplist')
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------       
    def test_file_line_selection(self):
        '''Test selecting lines from a file.'''
        f = open(TestItemUtil.__FILE, 'rb')
        num_lines = iu.ilen(f)
        assert_equal(num_lines, 3218, 'Wrong item count')
        
        f = open(TestItemUtil.__FILE, 'rb')
        a = list(iu.linerange(f, 0, num_lines, 100))
        assert_equal(len(a), 33, 'Wrong item range')
        assert_equal(a[0], 'SNP_A-2141031', 'Wrong item')
        assert_equal(a[-1], 'rs2015453', 'Wrong item')
        
    def test_segment_range(self):
        '''Test getting the start and end of equidistant intervals in an item collection.'''
        # Non-divisible case
        self.__test_segment_range(100, 7, 1)
        # Divisible case
        self.__test_segment_range(100, 10, 1)
        
        # Less trivial item list
        self.__test_segment_range(100, 7, 2)
        self.__test_segment_range(100, 10, 2)

    def test_groupby_with_range(self):
        '''Test the groupby_with_range() function.'''
        assert_equal(iu.groupby_with_range(''), [])
        assert_equal(iu.groupby_with_range('AAAABBBCCDAABBB'), 
                     [('A', (0, 3)),
                      ('B', (4, 6)),
                      ('C', (7, 8)),
                      ('D', (9, 9)),
                      ('A', (10, 11)),
                      ('B', (12, 14))])
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __test_segment_range(self, n, step, c):
        '''Test getting the start and end of equidistant intervals in an item collection.'''
        items = np.arange(0, n, c)
        segments = iu.segmentrange(items, step)
        assert_equal(segments[:, 0], range(0, n, c * step), 'Wrong start array')
        assert_equal(segments[:, 1],
                     range(c * (step - 1), n, c * step) + ([n - c] if np.mod(n, step) else []),
                     'Wrong end array')
