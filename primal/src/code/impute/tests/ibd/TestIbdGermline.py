'''
============================================================
Test the GERMLINE algorithm for finding IBD segments (in
particular in sibs). 

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

class TestIbdGermline(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------

    def setUp(self):
        '''Load single nuclear family test case.'''
        self.problem = im.io.read_npz(itu.SIB_FOUNDERS_STAGE3)
        self.family = self.problem.families(genotyped=False)[0]
        self.sibs = ig._filled_members(self.problem, self.family) 
        self.ibd_computer = ig.GermlineIbdComputer(PhaseParam())
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_match(self):
        '''Test hashing the rows of a haplotype matrix with the match() function.'''
        h_mat = ig._HapMatrix(self.problem, self.sibs)
        h = h_mat.h
        assert_equal(len(h_mat.snps), 3128, 'Wrong # filled snps')
        assert_equal(h.shape, (2 * len(self.sibs), len(h_mat.snps)),
                     'Wrong extracted hap matrix size')
        assert_equal(self.ibd_computer._match(h[:, 0:10]),
                     set([(1, 3), (3, 7), (0, 6), (0, 9), (1, 5),
                          (1, 7), (3, 5), (5, 7), (2, 4), (6, 9)]),
                     'Wrong match set')

    def test_ibs_segments(self):
        '''Test computing GERMLINE IBD segments.'''
        m = ig.ibd_germline(self.problem, self.sibs)
        assert_equal(m.length, 24, 'Wrong # IBD grouped segments')
        assert_segments_almost_equal(m,
                                     [((1   , 434), (17065079, 23614026, 6.549, 0), ((0, 1), (3, 1), (4, 1), (1, 1))),
                                      ((1   , 434), (17065079, 23614026, 6.549, 0), ((5, 1), (0, 0), (4, 0))),
                                      ((1   , 434), (17065079, 23614026, 6.549, 0), ((3, 0), (1, 0))),
                                      ((434 , 537), (23614026, 24732045, 1.118, 0), ((0, 1), (3, 1), (4, 1))),
                                      ((434 , 537), (23614026, 24732045, 1.118, 0), ((5, 1), (0, 0), (4, 0))),
                                      ((537 , 848), (24732045, 27321172, 2.589, 0), ((3, 0), (5, 1), (0, 0), (4, 0))),
                                      ((537 , 848), (24732045, 27321172, 2.589, 0), ((0, 1), (3, 1), (4, 1))),
                                      ((537 , 848), (24732045, 27321172, 2.589, 0), ((1, 1), (5, 0))),
                                      ((848 , 948), (27321172, 27671082, 0.350, 0), ((0, 1), (3, 1), (4, 1))),
                                      ((848 , 948), (27321172, 27671082, 0.350, 0), ((3, 0), (5, 1), (0, 0))),
                                      ((848 , 948), (27321172, 27671082, 0.350, 0), ((1, 1), (5, 0))),
                                      ((948 , 2481), (27671082, 44544028, 16.873, 0), ((0, 1), (3, 1), (4, 1))),
                                      ((948 , 2481), (27671082, 44544028, 16.873, 0), ((3, 0), (5, 1), (0, 0))),
                                      ((948 , 2481), (27671082, 44544028, 16.873, 0), ((1, 1), (5, 0))),
                                      ((948 , 2481), (27671082, 44544028, 16.873, 0), ((1, 0), (4, 0))),
                                      ((2481, 2583), (44544028, 45330235, 0.786, 0), ((0, 1), (3, 1), (4, 1))),
                                      ((2481, 2583), (44544028, 45330235, 0.786, 0), ((3, 0), (5, 1), (0, 0))),
                                      ((2481, 2583), (44544028, 45330235, 0.786, 0), ((1, 0), (4, 0))),
                                      ((2583, 3107), (45330235, 49669743, 4.340, 0), ((0, 1), (3, 1), (4, 1), (5, 0))),
                                      ((2583, 3107), (45330235, 49669743, 4.340, 0), ((3, 0), (5, 1), (0, 0))),
                                      ((2583, 3107), (45330235, 49669743, 4.340, 0), ((1, 0), (4, 0))),
                                      ((3107, 3218), (49669743, 51156933, 1.487, 0), ((0, 1), (3, 1), (4, 1), (1, 1), (5, 0))),
                                      ((3107, 3218), (49669743, 51156933, 1.487, 0), ((3, 0), (5, 1), (0, 0))),
                                      ((3107, 3218), (49669743, 51156933, 1.487, 0), ((1, 0), (4, 0)))],
                            full_data=True, decimal=3, err_msg='Wrong IBD segments')
        stats = np.array([(len(s.samples), s.length) for s in m])
        best_segment = np.lexsort((-stats[:, 1], -stats[:, 0]))[0]
        assert_equal(best_segment, 21, 'Wrong best segment (IBD set size + length)') 

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
