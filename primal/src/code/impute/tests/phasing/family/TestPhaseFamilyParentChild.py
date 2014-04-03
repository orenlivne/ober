'''
============================================================
Test phasing within nuclear families: Stage 2 (phased parent
-child IBD).

Created on October 25, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, impute as im
from numpy.ma.testutils import assert_equal
from impute.impute_test_util import assert_segments_almost_equal

class TestPhaseFamilyParentChild(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        # The way to load a pedigree in conjunction with a genotype set is to recode
        # its sample IDs to consecutive for easier access by phasers.
        self.phaser = im.phase_core.new_phaser_chain([im.phase_trivial.trivial_phaser(),
                                                      im.phase_family.family_phaser()])

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_global_ibd_segments(self):
        '''Debug error that appeared after stage 2 - child had wrong hap color at SNP 2960
        in the middle of a long IBD segment. This error doesn''t appear in this isolated test.
        It depends on the ordering of processing samples and segments.'''
        
        # Before stage 2: no IBD segments should exist
        problem = im.io.read_npz(im.itu.FAMILY4_STAGE1)
        im.itu.assert_problem_stats(problem, 38616, 35586, 22)
        ibd = problem.info.ibd
        assert_equal(ibd.length, 0, 'Before stage 2: no IBD segments should exist')
        assert_equal(problem.h[2955:2965, 1, 0], [2, 2, 2, 2, 2, 2, 2, 2, 2, 0], 'Wrong mother haplotype') 
        assert_equal(problem.h[2955:2965, 4, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 'Wrong child haplotype') 
        
        im.phase_family.family_phaser().run(problem, im.PhaseParam())
        
        # After: check IBD segment; (1, 0), (4, 1) should be IBD around SNP 2960
        assert_segments_almost_equal(ibd,
                                     [((0   , 176), (16484792, 18545634, 2.061, 0), ((2, 0), (0, 0))),
                                      ((187 , 2687), (18895227, 46242359, 27.347, 1), ((0, 1), (2, 0))),
                                      ((2687, 3216), (46336181, 51103692, 4.768, 0), ((2, 0), (0, 0))),
                                      ((0   , 3218), (16484792, 51156933, 34.672, 0), ((0, 1), (3, 0))),
                                      ((0   , 2470), (16484792, 44456692, 27.972, 0), ((0, 0), (4, 0))),
                                      ((2473, 3216), (44515806, 51103692, 6.588, 1), ((0, 1), (4, 0))),
                                      ((0   , 2985), (16484792, 48709188, 32.224, 0), ((0, 1), (5, 0))),
                                      ((2993, 3216), (48742097, 51103692, 2.362, 0), ((0, 0), (5, 0))),
                                      ((0   , 3218), (16484792, 51156933, 34.672, 0), ((1, 1), (2, 1))),
                                      ((4   , 36), (17087656, 17434521, 0.347, 0), ((3, 1), (1, 1))),
                                      ((42  , 2035), (17587680, 37662436, 20.075, 1), ((1, 0), (3, 1))),
                                      ((2047, 3217), (37902926, 51140316, 13.237, 0), ((3, 1), (1, 1))),
                                      ((11  , 804), (17285049, 27091750, 9.807, 0), ((4, 1), (1, 1))),
                                      ((823 , 3217), (27200942, 51140316, 23.939, 0), ((1, 0), (4, 1))),
                                      ((11  , 14), (17285049, 17307742, 0.023, 0), ((5, 1), (1, 0))),
                                      ((31  , 3217), (17415572, 51140316, 33.725, 1), ((5, 1), (1, 1)))],
                                     decimal=3, err_msg='Wrong IBD segments') 
        im.itu.assert_problem_stats(problem, 38616, 38576, 30)
        assert_equal(problem.h[2955:2965, 1, 0], [2, 2, 2, 2, 2, 2, 2, 2, 2, 1], 'Wrong mother haplotype') 
        assert_equal(problem.h[2955:2965, 4, 1], [2, 2, 2, 2, 2, 2, 2, 2, 2, 1], 'Wrong child haplotype') 

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
