'''
============================================================
Test new haplotype coloring algorithm (matrix-based). 

Created on January 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, impute as im, itertools as it
from numpy.ma.testutils import assert_equal

class TestHapColoring(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        pass
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_quasi_founder_sibs(self):
        '''Check segment coloring for a quasi-founder family with 5 children. Only the 4
        sibs that were POO-aligned are considered.'''
        sibs = im.itu.QF_FAMILY5_SIBS4
        segments = im.segment.SegmentSet.load(open(im.itu.QF_FAMILY5_SIBS4_SEGMENTS, 'rb'))
        
        pa = im.color.hap_color.hap_colors(list(it.product(sorted(sibs), im.constants.ALLELES)), segments, max_colors=4)
        assert_equal(pa.haps, [(250, 0), (250, 1), (597, 0), (597, 1), (631, 0), (631, 1), (1028, 0), (1028, 1)], 'Wrong palette haplotype list')
        assert_equal(pa.regions, [(0, 15), (15, 59), (59, 104), (104, 518), (518, 519), (519, 635), (635, 748), (748, 759), (759, 968), (968, 969), (969, 3075), (3075, 3108), (3108, 3217)], 'Wrong palette region list')
        assert_equal(pa.color, [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0],
                                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [2, 2, 0, 0, 0, 0, 2, 1, 1, 2, 2, 2, 2],
                                [2, 2, 0, 0, 0, 0, 2, 2, 2, 4, 1, 1, 1],
                                [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [2, 2, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]], 'Wrong palette colors')
        assert_equal(pa.color_sequence(), [[(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9), (0, 10), (2, 11), (0, 12)], [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (1, 11), (1, 12)], [(3, 0), (3, 1), (0, 2), (0, 3), (0, 4), (0, 5), (3, 6), (4, 7), (4, 8), (3, 9), (3, 10), (3, 11), (3, 12)], [(5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (5, 10), (5, 11), (0, 12)]], 'Wrong paternal hap color sequence')

    def test_non_founder_sibs(self):
        '''Check segment coloring for a non-founder family with 4 children.'''
        sibs = im.itu.NF_FAMILY_SIBS
        segments = im.segment.SegmentSet.load(open(im.itu.NF_FAMILY_SEGMENTS, 'rb'))
        
        pa = im.color.hap_color.hap_colors(list(it.product(sorted(sibs), im.constants.ALLELES)), segments, max_colors=4)
        assert_equal(pa.color, [[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                                [0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 3, 3, 3, 3, 3, 3],
                                [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 4, 4, 2],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 2, 2, 4, 1],
                                [0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], 'Wrong palette colors')

#---------------------------------------------
# Private Methods
#---------------------------------------------
