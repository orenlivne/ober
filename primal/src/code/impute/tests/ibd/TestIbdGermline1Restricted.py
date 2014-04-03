'''
============================================================
Test the GERMLINE algorithm for finding IBD segments among
sibs with one genotyped parent.

Created on September 14, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, impute as im
from numpy.ma.testutils import assert_equal, assert_almost_equal
from impute import impute_test_util as itu
from impute.ibd import ibd_germline as ig
from impute.tools.param import PhaseParam
from impute.impute_test_util import assert_segments_almost_equal
from impute.data.constants import MATERNAL
from hashable import Hashable

class TestIbdGermline1Restricted(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------

    def setUp(self):
        '''Load single nuclear family test case.'''
        # Remove a key child to make problem more interesting for the IBD algorithm
        self.problem = im.io.read_npz(itu.FAMILY945_ONE_PARENT_STAGE2).remove_nodes([2])
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
        
    def test_match_both_parents(self):
        '''Test hashing the rows of a haplotype matrix with the match() function. Here H is n x s.'''
        h_mat = ig._HapMatrix(self.problem, self.sibs, parent_type=MATERNAL)
        h = h_mat.h
        assert_equal(len(h_mat.snps), 2929, 'Wrong # filled snps')
        assert_equal(h.shape, (len(self.sibs), len(h_mat.snps)), 'Wrong extracted hap matrix size')
        #h_mat.print_info()
        assert_equal(self.ibd_computer._match(h[:, 0:10]), set([(1, 2)]), 'Wrong match set')
        assert_equal(self.ibd_computer._match(h[:, 2100:2400]), set([(0, 1), (1, 2), (0, 2)]), 'Wrong match set')

    def test_ibd_segments_both_parents(self):
        '''Test computing GERMLINE IBD segments.'''
        h_mat = ig._HapMatrix(self.problem, self.sibs)
        m = self.ibd_computer.ibd_segments(h_mat)
        assert_segments_almost_equal(m,
                                     [((0   , 344), (16484792, 21449028, 4.964, 0), ((2, 0), (4, 0))),
                                      ((0   , 344), (16484792, 21449028, 4.964, 0), ((1, 0), (4, 0))),
                                      ((0   , 780), (16484792, 26920270, 10.435, 0), ((4, 1), (2, 1))),
                                      ((2145, 2996), (39608193, 48759228, 9.151, 0), ((1, 1), (2, 1))),
                                      ((2145, 2996), (39608193, 48759228, 9.151, 0), ((4, 1), (2, 1))),
                                      ((885 , 3218), (27425790, 51156933, 23.731, 0), ((4, 1), (1, 1))),
                                      ((0   , 3218), (16484792, 51156933, 34.672, 0), ((2, 0), (1, 0)))],
                                     full_data=True, decimal=3, err_msg='Wrong IBD segments, raw')
        # A transitive-logic test test to see that we don't miss any intervals with GERMLINE
        m.group_to_disjoint(False) 
        assert_segments_almost_equal(m,
                                     [((0   , 344), (16484792, 21449028, 4.964, 0), ((2, 0), (1, 0), (4, 0))),
                                      ((0   , 344), (16484792, 21449028, 4.964, 0), ((4, 1), (2, 1))),
                                      ((344 , 780), (21449028, 26920270, 5.471, 0), ((2, 0), (1, 0))),
                                      ((344 , 780), (21449028, 26920270, 5.471, 0), ((4, 1), (2, 1))),
                                      ((780 , 885), (26920270, 27425790, 0.506, 0), ((2, 0), (1, 0))),
                                      ((885 , 2145), (27425790, 39608193, 12.182, 0), ((2, 0), (1, 0))),
                                      ((885 , 2145), (27425790, 39608193, 12.182, 0), ((4, 1), (1, 1))),
                                      ((2145, 2996), (39608193, 48759228, 9.151, 0), ((4, 1), (1, 1), (2, 1))),
                                      ((2145, 2996), (39608193, 48759228, 9.151, 0), ((2, 0), (1, 0))),
                                      ((2996, 3218), (48759228, 51156933, 2.398, 0), ((2, 0), (1, 0))),
                                      ((2996, 3218), (48759228, 51156933, 2.398, 0), ((4, 1), (1, 1)))],
                                     full_data=True, decimal=3, err_msg='Wrong IBD segments, grouped')

        stats = np.array([(len(s.samples), s.length) for s in m])
        best_segment = np.lexsort((-stats[:, 1], -stats[:, 0]))[0]
        assert_equal(best_segment, 7, 'Wrong best segment (IBD set size + length)') 
        assert_almost_equal(m[best_segment].length, 9.15, decimal=2, err_msg='Wrong best segment (IBD set size + length)')
        assert_equal(m[best_segment].samples, set([(4, 1), (1, 1), (2, 1)]), err_msg='Wrong best segment (IBD set size + length)')

    def test_match_maternal(self):
        '''Test hashing the rows of a haplotype matrix with the match() function. Here H is n x s.'''
        h_mat = ig._HapMatrix(self.problem, self.sibs, parent_type=MATERNAL)
        h = h_mat.h
        assert_equal(len(h_mat.snps), 2929, 'Wrong # filled snps')
        assert_equal(h.shape, (len(self.sibs), len(h_mat.snps)), 'Wrong extracted hap matrix size')
        assert_equal(self.ibd_computer._match(h[:, 0:10]), set([(1, 2)]), 'Wrong match set')

    def test_ibd_segments_maternal(self):
        '''Test computing GERMLINE IBD segments among maternal haps only.'''
        pass
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
