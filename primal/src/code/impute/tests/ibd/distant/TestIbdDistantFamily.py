'''
============================================================
Test finding IBD segments within a nuclear families using
the general IBD-posterior strategy. 

Created on July 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, impute as im
from impute.tools.param import PhaseParam
from impute.ibd.distant.ibd_distant import ibd_segments_with_relatives
from numpy.ma.testutils import assert_equal
from impute.tools.genotype_tools import genotyped_children

class TestIbdDistantFamily(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        
        # Load test data ready from previous phasing stagees
        self.problem = im.io.read_npz(im.itu.FAMILY_TOO_ZEROED_STAGE4) 
        self.family = self.problem.first_family

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_ibd_segments_sib_pair(self):
        '''Test calculating distant IBD segments between a pair of sibs;
        compare with IBD segments based on nucelar family info'''
        segment_set = ibd_segments_with_relatives(self.problem, 3, 2,
                                                  PhaseParam(id_coef_file=im.itu.ID_COEF_FILE),
                                                  im.ibd_hmm.prob_ibd_hmm)

        im.itu.assert_segments_almost_equal(segment_set,
                                            [((1412, 3218), (32992389, 51156934, 18.165, 1), ((3, 0), (2, 0))),
                                             ((241 , 451), (19643555, 23817486, 4.174, 1), ((2, 0), (3, 1))),
                                             ((0   , 600), (16484792, 25444874, 8.960, 1), ((3, 1), (2, 1))),
                                             ((2650, 3218), (45892433, 51156934, 5.265, 1), ((3, 1), (2, 1)))],
                                             full_data=True, decimal=3, err_msg='Wrong IBD segments')

    def test_ibd_segments_sib_pair2(self):
        '''Test calculating distant IBD segments against a single surrogate parent;
        compare with IBD segments based on nucelar family info.'''
        segment_set = ibd_segments_with_relatives(self.problem, 3, 5,
                                                  PhaseParam(id_coef_file=im.itu.ID_COEF_FILE, max_path_length=2),
                                                  im.ibd_hmm.prob_ibd_hmm)

        im.itu.assert_segments_almost_equal(segment_set,
                                            [((1278, 1337), (31217345, 32331594, 1.114, 1), ((3, 0), (5, 0))),
                                             ((2276, 2363), (42344297, 43527681, 1.183, 1), ((3, 1), (5, 0))),
                                             ((3138, 3206), (49803008, 50837224, 1.034, 1), ((3, 1), (5, 0))),
                                             ((603 , 3218), (25453554, 51156934, 25.703, 1), ((5, 1), (3, 1)))],
                                             full_data=True, decimal=3, err_msg='Wrong IBD segments')

    def test_ibd_parent_vs_all_children(self):
        '''Test calculating distant IBD segments against all surrogate parents;
        compare with IBD segments based on nucelar family info.'''
#        segment_set = ibd_segments_with_surrogate_parents(self.problem, 0,
#                                                          PhaseParam(margin=0., surrogate_parent_fill_threshold=0.9,
#                                                                     max_path_length=2, debug=True),
#                                                          prob_ibd_calculator=im.ibd_hmm.prob_ibd_hmm,
#                                                          is_i_phased=True)
        # Turn off kinship-based POO determination IBD segment computation since we don't have
        # a complete pedigree here 
        segment_set = ibd_segments_with_relatives(self.problem, 0,
                                                  genotyped_children(self.problem, self.problem.first_family),
                                                  PhaseParam(id_coef_file=im.itu.ID_COEF_FILE, max_path_length=2),
                                                  im.ibd_hmm.prob_ibd_hmm, use_kinship=False)
        im.itu.assert_segments_almost_equal(segment_set,
                                            [((0   , 1420), (16484792, 33032458, 16.548, 1), ((2, 0), (0, 0))),
                                             ((241 , 446), (19643555, 23761236, 4.118, 1), ((0, 0), (2, 1))),
                                             ((2215, 2270), (40876234, 42241372, 1.365, 1), ((0, 0), (2, 1))),
                                             ((3138, 3206), (49803008, 50837224, 1.034, 1), ((0, 0), (2, 1))),
                                             ((1278, 1337), (31217345, 32331594, 1.114, 1), ((0, 1), (2, 0))),
                                             ((1411, 3218), (32978753, 51156934, 18.178, 1), ((0, 1), (2, 0))),
                                             ((1278, 1337), (31217345, 32331594, 1.114, 1), ((3, 0), (0, 0))),
                                             ((241 , 451), (19643555, 23817486, 4.174, 1), ((0, 0), (3, 1))),
                                             ((2276, 2363), (42344297, 43527681, 1.183, 1), ((0, 0), (3, 1))),
                                             ((3138, 3206), (49803008, 50837224, 1.034, 1), ((0, 0), (3, 1))),
                                             ((0   , 3218), (16484792, 51156934, 34.672, 1), ((0, 1), (3, 0))),
                                             ((1278, 1337), (31217345, 32331594, 1.114, 1), ((0, 0), (4, 0))),
                                             ((2552, 3218), (45011952, 51156934, 6.145, 1), ((0, 0), (4, 0))),
                                             ((385 , 451), (22583252, 23817486, 1.234, 1), ((0, 0), (4, 1))),
                                             ((2215, 2270), (40876234, 42241372, 1.365, 1), ((0, 0), (4, 1))),
                                             ((0   , 2571), (16484792, 45231758, 28.747, 1), ((0, 1), (4, 0))),
                                             ((0   , 3218), (16484792, 51156934, 34.672, 1), ((0, 0), (5, 0))),
                                             ((385 , 451), (22583252, 23817486, 1.234, 1), ((5, 1), (0, 0))),
                                             ((2276, 2363), (42344297, 43527681, 1.183, 1), ((5, 1), (0, 0))),
                                             ((3138, 3206), (49803008, 50837224, 1.034, 1), ((5, 1), (0, 0))),
                                             ((1278, 1337), (31217345, 32331594, 1.114, 1), ((0, 1), (5, 0))),
                                             ((0   , 805), (16484792, 27119061, 10.634, 1), ((0, 0), (6, 0))),
                                             ((1278, 1337), (31217345, 32331594, 1.114, 1), ((0, 0), (6, 0))),
                                             ((241 , 347), (19643555, 22015144, 2.372, 1), ((6, 1), (0, 0))),
                                             ((385 , 451), (22583252, 23817486, 1.234, 1), ((6, 1), (0, 0))),
                                             ((2276, 2363), (42344297, 43527681, 1.183, 1), ((6, 1), (0, 0))),
                                             ((3138, 3206), (49803008, 50837224, 1.034, 1), ((6, 1), (0, 0))),
                                             ((762 , 3218), (26836780, 51156934, 24.320, 1), ((0, 1), (6, 0)))],
                                             full_data=True, decimal=3, err_msg='Wrong IBD segments')

    def test_ibd_from_nuclear_family(self):
        '''Test getting the correct IBD picture using nuclear family phasing.'''
        p = im.io.read_npz(im.itu.FAMILY_TOO_ZEROED_STAGE2) 
        ibd = p.info.ibd
        assert_equal(ibd.length, 17, 'Unexpected final # of IBD segments before grouping')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
