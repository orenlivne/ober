'''
============================================================
Test finding IBD segments between haplotypes of nuclear
family members. Uses HMM-Hap. 

Created on January 27, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, impute as im
# from numpy.ma.testutils import assert_equal

class TestIbdDistantHapFamily(unittest.TestCase):
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
        self.problem = im.io.read_npz(im.itu.FAMILY_TOO_ZEROED_STAGE5) 
        self.family = self.problem.first_family

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_ibd_segments_sib_pair(self):
        '''Test calculating distant IBD segments between every hap pair within a pair of phased sibs.
        Compare with IBD segments based on nucelar family info.'''
        # Segments should contain segments obtained by the genotype-genotye-HMM to validate phasing:
        #        [((1412, 3218), (32992389, 51156934, 18.165, 1), ((3, 0), (2, 0))),
        #         ((241 , 451), (19643555, 23817486, 4.174, 1), ((2, 0), (3, 1))),
        #         ((0   , 454), (16484792, 23834889, 7.350, 1), ((3, 1), (2, 1))),
        #         ((2650, 3218), (45892433, 51156934, 5.265, 1), ((3, 1), (2, 1)))],
        expected_segments = [((1412, 3217), (32992389, 51156933, 18.165, 0), ((3, 0), (2, 0))),
                             ((238 , 453), (19581946, 23826675, 4.245, 0), ((2, 0), (3, 1))),
                             ((0   , 629), (16484792, 25608548, 9.124, 0), ((3, 1), (2, 1))),
                             ((2661, 3217), (45972017, 51156933, 5.185, 0), ((3, 1), (2, 1)))]
        # Serial test
        segment_set = im.ih.between_samples_segments(self.problem, [3], [2], im.PhaseParam(kinship_file=im.itu.KINSHIP_FILE, debug=False))
        im.itu.assert_segments_almost_equal(segment_set, expected_segments, full_data=True, decimal=3, err_msg='Wrong IBD segments')

        # Parallel test
        segment_set = im.ih.between_samples_segments(self.problem, [3], [2], im.PhaseParam(kinship_file=im.itu.KINSHIP_FILE, debug=False), num_processes=3)
        im.itu.assert_segments_almost_equal(segment_set, expected_segments, full_data=True, decimal=3, err_msg='Wrong IBD segments')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
