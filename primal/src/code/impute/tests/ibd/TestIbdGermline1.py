'''
============================================================
Test the GERMLINE algorithm for finding IBD segments among
sibs with one genotyped parent.

Created on September 14, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, impute as im
from numpy.ma.testutils import assert_equal
from impute import impute_test_util as itu
from impute.ibd import ibd_germline as ig
from impute.tools.param import PhaseParam
from impute.impute_test_util import assert_segments_almost_equal
from impute.data.constants import MATERNAL
from hashable import Hashable

class TestIbdGermline1(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------

    def setUp(self):
        '''Load single nuclear family test case.'''
        self.problem = im.io.read_npz(itu.FAMILY945_ONE_PARENT_STAGE2)
        self.family = self.problem.families(genotyped=False)[0]
        self.sibs = ig._filled_members(self.problem, self.family) 
        self.ibd_computer = ig.GermlineIbdComputer(PhaseParam())
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_hash(self):
        '''Test hashing the rows of a haplotype matrix with the match() function. Here H is n x s.'''
        h_mat = ig._HapMatrix(self.problem, self.sibs, parent_type=MATERNAL)
        h = h_mat.h
        h = h[:, 0:10]
        assert_equal(h[1], h[2], 'These two test rows should be equal')
        def hash_array(x): x.tostring()
        assert_equal(hash_array(h[1]), hash_array(h[2]), 'Test objects should be equal, if not, rerun test')

        d = {}
        for (i, hi) in enumerate(h):
            d.setdefault(Hashable(hi), []).append(i)
        assert_equal(len(d), 2, 'Hashing haps into dictionary did not work')
        
    def test_match(self):
        '''Test hashing the rows of a haplotype matrix with the match() function. Here H is n x s.'''
        h_mat = ig._HapMatrix(self.problem, self.sibs, parent_type=MATERNAL)
        h = h_mat.h
        assert_equal(len(h_mat.snps), 2915, 'Wrong # filled snps')
        assert_equal(h.shape, (len(self.sibs), len(h_mat.snps)),
                     'Wrong extracted hap matrix size')
        assert_equal(self.ibd_computer._match(h[:, 0:10]),
                     set([(1, 2), (1, 3), (2, 3)]),
                     'Wrong match set')

    def test_ibd_segments(self):
        '''Test computing GERMLINE IBD segments.'''
        h_mat = ig._HapMatrix(self.problem, self.sibs)
        m = self.ibd_computer.ibd_segments(h_mat)
        m.group_to_disjoint()
        assert_segments_almost_equal(m,
                                    [((0   , 235), (16484792, 19567818, 3.083, 0), ((3, 0), (1, 0), (5, 0))),
                                     ((0   , 235), (16484792, 19567818, 3.083, 0), ((5, 1), (3, 1), (2, 1))),
                                     ((235 , 345), (19567818, 21809185, 2.241, 0), ((3, 0), (1, 0), (5, 0))),
                                     ((235 , 345), (19567818, 21809185, 2.241, 0), ((5, 1), (3, 1))),
                                     ((345 , 450), (21809185, 23817240, 2.008, 0), ((3, 0), (1, 0))),
                                     ((345 , 450), (21809185, 23817240, 2.008, 0), ((5, 1), (3, 1))),
                                     ((345 , 450), (21809185, 23817240, 2.008, 0), ((1, 1), (2, 1))),
                                     ((450 , 781), (23817240, 26924312, 3.107, 0), ((3, 0), (1, 0))),
                                     ((450 , 781), (23817240, 26924312, 3.107, 0), ((5, 1), (3, 1))),
                                     ((450 , 781), (23817240, 26924312, 3.107, 0), ((1, 1), (2, 1))),
                                     ((450 , 781), (23817240, 26924312, 3.107, 0), ((2, 0), (5, 0))),
                                     ((781 , 886), (26924312, 27427698, 0.503, 0), ((3, 0), (1, 0))),
                                     ((781 , 886), (26924312, 27427698, 0.503, 0), ((1, 1), (2, 1))),
                                     ((781 , 886), (26924312, 27427698, 0.503, 0), ((2, 0), (5, 0))),
                                     ((886 , 2157), (27427698, 39797178, 12.369, 0), ((5, 1), (1, 1), (2, 1))),
                                     ((886 , 2157), (27427698, 39797178, 12.369, 0), ((3, 0), (1, 0))),
                                     ((886 , 2157), (27427698, 39797178, 12.369, 0), ((2, 0), (5, 0))),
                                     ((2157, 2791), (39797178, 47440321, 7.643, 0), ((5, 1), (3, 1), (1, 1), (2, 1))),
                                     ((2157, 2791), (39797178, 47440321, 7.643, 0), ((3, 0), (1, 0))),
                                     ((2157, 2791), (39797178, 47440321, 7.643, 0), ((2, 0), (5, 0))),
                                     ((2791, 2901), (47440321, 48282594, 0.842, 0), ((5, 1), (3, 1), (1, 1), (2, 1))),
                                     ((2791, 2901), (47440321, 48282594, 0.842, 0), ((3, 0), (1, 0))),
                                     ((2901, 3012), (48282594, 48855412, 0.573, 0), ((5, 1), (3, 1), (1, 1), (2, 1))),
                                     ((2901, 3012), (48282594, 48855412, 0.573, 0), ((3, 0), (2, 0), (1, 0))),
                                     ((3012, 3218), (48855412, 51156933, 2.302, 0), ((3, 0), (2, 0), (1, 0))),
                                     ((3012, 3218), (48855412, 51156933, 2.302, 0), ((5, 1), (1, 1), (2, 1)))],
                            full_data=True, decimal=3, err_msg='Wrong IBD segments')
        stats = np.array([(len(s.samples), s.length) for s in m])
        best_segment = np.lexsort((-stats[:, 1], -stats[:, 0]))[0]
        #print np.lexsort((-stats[:, 1], -stats[:, 0]))
        # print m[best_segment]
        assert_equal(best_segment, 17, 'Wrong best segment (IBD set size + length)') 

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
