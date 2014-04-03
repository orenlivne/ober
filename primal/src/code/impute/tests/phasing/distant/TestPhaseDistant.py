'''
============================================================
Test phasing within nuclear families.

Created on July 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, impute as im
from impute.impute_test_util import assert_segments_almost_equal
from numpy.ma.testutils import assert_almost_equal, assert_equal

class TestPhaseDistant(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        self.problem = im.io.read_npz(im.itu.NBHRS1298_STAGE4_WITH_PARENTS)

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_families(self):
        '''Test that nuclear families computation does not crash for a disconnected pedigree
        such as NBHRS1298_STAGE4.'''
        f = list(self.problem.families())
        assert_equal(len(f), 2, 'Wrong number of nuclear families')
        # self.assertFalse(f, 'No nuclear families should be identified') # there are some now since we added parents to the nbhr set

    def test_phasing_distant(self):
        '''Check phasing a single sample based on distant relatives.'''
        proband = 14  # Unphased sample index
        phaser = im.phase_core.new_phaser_chain([im.phase_distant.distant_phaser(single_sample=proband)])
        im.itu.assert_size_equals(self.problem.genotype, 3218, 28)
        h = self.problem.haplotype
        im.itu.assert_problem_stats(self.problem, 180208, 150235, 114)
        assert_almost_equal(h.fill_fraction(sample=proband), 0.46, 2, 'Unexpected pre-phasing fill %')
        phaser.run(self.problem, im.PhaseParam(target_fill=0.95, max_path_length=7))
        im.itu.assert_problem_stats(self.problem, 180208, 153675, 114)
        assert_almost_equal(h.fill_fraction(sample=proband), 0.943, 2, 'Unexpected post-phasing fill %')
        
    def test_ibd_segments(self):
        '''Test locating IBD segments between distant individuals..'''
        segment_set = im.ibd_distant.ibd_segments_with_surrogate_parents(self.problem, 14, 1, 6, params=im.PhaseParam(debug=False),
                                                                         prob_ibd_calculator=im.ibd_hmm.prob_ibd_hmm)
        segment_set.group_to_disjoint()
        assert_segments_almost_equal(segment_set,
[((0   , 76), (16484792, 17844295, 1.360, 0), ((16, 1), (16, 0), (23, 0), (14, 0), (23, 1), (1, 0), (1, 1))),
 ((0   , 76), (16484792, 17844295, 1.360, 0), ((5, 1), (14, 1), (5, 0))),
 ((76  , 90), (17844295, 17909524, 0.065, 0), ((16, 1), (23, 0), (14, 0), (23, 1), (1, 0), (1, 1))),
 ((76  , 90), (17844295, 17909524, 0.065, 0), ((5, 1), (14, 1), (5, 0))),
 ((90  , 95), (17909524, 17930263, 0.021, 0), ((16, 1), (16, 0), (23, 0), (14, 0), (23, 1), (1, 0), (1, 1))),
 ((90  , 95), (17909524, 17930263, 0.021, 0), ((5, 1), (14, 1), (5, 0))),
 ((95  , 96), (17930263, 17948473, 0.018, 0), ((16, 1), (16, 0), (14, 0), (23, 1), (1, 0), (1, 1))),
 ((95  , 96), (17930263, 17948473, 0.018, 0), ((5, 1), (14, 1), (5, 0))),
 ((96  , 100), (17948473, 17987962, 0.039, 0), ((16, 1), (16, 0), (14, 0), (1, 1))),
 ((100 , 226), (17987962, 19486391, 1.498, 0), ((16, 1), (16, 0), (23, 0), (14, 0), (23, 1), (1, 0), (1, 1))),
 ((100 , 226), (17987962, 19486391, 1.498, 0), ((5, 1), (14, 1), (5, 0))),
 ((226 , 251), (19486391, 19787736, 0.301, 0), ((16, 0), (23, 0), (14, 0), (23, 1), (1, 0), (1, 1))),
 ((226 , 251), (19486391, 19787736, 0.301, 0), ((5, 1), (14, 1), (5, 0))),
 ((251 , 275), (19787736, 20075858, 0.288, 0), ((16, 1), (16, 0), (23, 0), (14, 0), (23, 1), (1, 0), (1, 1))),
 ((251 , 275), (19787736, 20075858, 0.288, 0), ((5, 1), (14, 1), (5, 0))),
 ((275 , 279), (20075858, 20742450, 0.667, 0), ((16, 1), (1, 0), (14, 0))),
 ((279 , 300), (20742450, 20993519, 0.251, 0), ((16, 1), (23, 0), (14, 0), (23, 1), (1, 0), (1, 1))),
 ((279 , 300), (20742450, 20993519, 0.251, 0), ((5, 1), (14, 1), (5, 0))),
 ((300 , 380), (20993519, 22554306, 1.561, 0), ((16, 1), (16, 0), (23, 0), (14, 0), (23, 1), (1, 0), (1, 1))),
 ((300 , 380), (20993519, 22554306, 1.561, 0), ((5, 1), (14, 1), (5, 0))),
 ((380 , 383), (22554306, 22555078, 0.001, 0), ((16, 1), (1, 0), (23, 1), (1, 1), (14, 0))),
 ((380 , 383), (22554306, 22555078, 0.001, 0), ((14, 1), (5, 0))),
 ((383 , 385), (22555078, 22583252, 0.028, 0), ((16, 1), (16, 0), (14, 0), (23, 1), (1, 0), (1, 1))),
 ((383 , 385), (22555078, 22583252, 0.028, 0), ((14, 1), (5, 0))),
 ((385 , 401), (22583252, 23275227, 0.692, 0), ((16, 1), (16, 0), (23, 0), (14, 0), (23, 1), (1, 0), (1, 1))),
 ((385 , 401), (22583252, 23275227, 0.692, 0), ((5, 1), (14, 1), (5, 0))),
 ((401 , 444), (23275227, 23753059, 0.478, 0), ((0, 1), (16, 1), (16, 0), (23, 0), (14, 0), (23, 1), (24, 1), (1, 0), (1, 1))),
 ((401 , 444), (23275227, 23753059, 0.478, 0), ((5, 1), (14, 1), (5, 0))),
 ((444 , 454), (23753059, 23834889, 0.082, 0), ((0, 1), (16, 1), (23, 0), (14, 0), (23, 1), (24, 1), (1, 0), (1, 1))),
 ((444 , 454), (23753059, 23834889, 0.082, 0), ((5, 1), (14, 1), (5, 0))),
 ((454 , 466), (23834889, 23877733, 0.043, 0), ((0, 1), (23, 0), (14, 0), (23, 1), (24, 1), (1, 0), (1, 1))),
 ((454 , 466), (23834889, 23877733, 0.043, 0), ((5, 1), (14, 1), (5, 0))),
 ((466 , 499), (23877733, 24209759, 0.332, 0), ((0, 1), (27, 1), (0, 0), (23, 0), (14, 0), (23, 1), (24, 1), (1, 0), (24, 0), (1, 1))),
 ((466 , 499), (23877733, 24209759, 0.332, 0), ((5, 1), (14, 1), (5, 0))),
 ((499 , 570), (24209759, 25205719, 0.996, 0), ((0, 1), (27, 1), (0, 0), (14, 0), (23, 1), (24, 1), (1, 0), (24, 0), (1, 1))),
 ((499 , 570), (24209759, 25205719, 0.996, 0), ((14, 1), (5, 0))),
 ((570 , 582), (25205719, 25305370, 0.100, 0), ((0, 1), (27, 0), (27, 1), (0, 0), (14, 0), (23, 1), (24, 1), (1, 0), (24, 0), (1, 1))),
 ((570 , 582), (25205719, 25305370, 0.100, 0), ((14, 1), (5, 0))),
 ((582 , 584), (25305370, 25319244, 0.014, 0), ((0, 1), (27, 0), (27, 1), (0, 0), (14, 0), (23, 1), (24, 1), (24, 0))),
 ((582 , 584), (25305370, 25319244, 0.014, 0), ((14, 1), (5, 0))),
 ((584 , 674), (25319244, 26234080, 0.915, 0), ((0, 1), (27, 0), (27, 1), (7, 1), (0, 0), (14, 0), (23, 1), (24, 1), (24, 0))),
 ((584 , 674), (25319244, 26234080, 0.915, 0), ((14, 1), (5, 0))),
 ((674 , 679), (26234080, 26250271, 0.016, 0), ((0, 1), (27, 1), (7, 1), (0, 0), (14, 0), (23, 1), (24, 1), (24, 0))),
 ((674 , 679), (26234080, 26250271, 0.016, 0), ((14, 1), (5, 0))),
 ((679 , 700), (26250271, 26380870, 0.131, 0), ((0, 1), (0, 0), (26, 1), (7, 1), (27, 1), (14, 0), (23, 1), (24, 1), (24, 0))),
 ((679 , 700), (26250271, 26380870, 0.131, 0), ((14, 1), (5, 0))),
 ((700 , 710), (26380870, 26496622, 0.116, 0), ((0, 1), (14, 0), (26, 1), (7, 1), (24, 1), (23, 1))),
 ((700 , 710), (26380870, 26496622, 0.116, 0), ((14, 1), (5, 0))),
 ((710 , 1008), (26496622, 28019908, 1.523, 0), ((0, 1), (9, 1), (26, 1), (7, 1), (24, 1), (23, 1), (14, 0))),
 ((710 , 1008), (26496622, 28019908, 1.523, 0), ((14, 1), (5, 0))),
 ((1008, 1032), (28019908, 28120707, 0.101, 0), ((0, 1), (7, 0), (9, 1), (26, 1), (7, 1), (14, 0), (23, 1), (24, 1))),
 ((1008, 1032), (28019908, 28120707, 0.101, 0), ((14, 1), (5, 0))),
 ((1032, 1039), (28120707, 28169949, 0.049, 0), ((0, 1), (7, 0), (9, 1), (26, 1), (7, 1), (23, 0), (14, 0), (23, 1), (24, 1))),
 ((1032, 1039), (28120707, 28169949, 0.049, 0), ((5, 1), (14, 1), (5, 0))),
 ((1039, 1088), (28169949, 29130300, 0.960, 0), ((0, 1), (9, 0), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (23, 0), (14, 0), (23, 1), (24, 1))),
 ((1039, 1088), (28169949, 29130300, 0.960, 0), ((5, 1), (14, 1), (5, 0))),
 ((1088, 1091), (29130300, 29145539, 0.015, 0), ((0, 1), (9, 0), (0, 0), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (23, 0), (14, 0), (23, 1), (24, 1))),
 ((1088, 1091), (29130300, 29145539, 0.015, 0), ((5, 1), (14, 1), (5, 0))),
 ((1091, 1123), (29145539, 29465677, 0.320, 0), ((0, 1), (9, 0), (0, 0), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (14, 0), (23, 1), (24, 1))),
 ((1091, 1123), (29145539, 29465677, 0.320, 0), ((14, 1), (5, 0))),
 ((1123, 1146), (29465677, 29661341, 0.196, 0), ((0, 1), (9, 0), (0, 0), (9, 1), (26, 1), (7, 1), (26, 0), (14, 0), (23, 1), (24, 1))),
 ((1123, 1146), (29465677, 29661341, 0.196, 0), ((14, 1), (5, 0))),
 ((1146, 1160), (29661341, 29834766, 0.173, 0), ((0, 1), (9, 0), (0, 0), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (14, 0), (23, 1), (24, 1))),
 ((1146, 1160), (29661341, 29834766, 0.173, 0), ((14, 1), (5, 0))),
 ((1160, 1207), (29834766, 30433373, 0.599, 0), ((0, 1), (9, 0), (0, 0), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (14, 0), (23, 1), (24, 1), (1, 1))),
 ((1160, 1207), (29834766, 30433373, 0.599, 0), ((14, 1), (5, 0))),
 ((1207, 1216), (30433373, 30575812, 0.142, 0), ((0, 1), (9, 0), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (14, 0), (23, 1), (24, 1), (1, 1))),
 ((1207, 1216), (30433373, 30575812, 0.142, 0), ((14, 1), (5, 0))),
 ((1216, 1283), (30575812, 31288769, 0.713, 0), ((0, 1), (9, 0), (0, 0), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (14, 0), (23, 1), (24, 1), (1, 1))),
 ((1216, 1283), (30575812, 31288769, 0.713, 0), ((14, 1), (5, 0))),
 ((1283, 1284), (31288769, 31323900, 0.035, 0), ((0, 1), (0, 0), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (14, 0), (23, 1), (24, 1), (1, 1))),
 ((1283, 1284), (31288769, 31323900, 0.035, 0), ((14, 1), (5, 0))),
 ((1284, 1298), (31323900, 31546951, 0.223, 0), ((0, 1), (9, 0), (0, 0), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (14, 0), (23, 1), (24, 1), (1, 1))),
 ((1284, 1298), (31323900, 31546951, 0.223, 0), ((14, 1), (5, 0))),
 ((1298, 1326), (31546951, 32210300, 0.663, 0), ((0, 1), (9, 0), (0, 0), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (23, 0), (14, 0), (23, 1), (24, 1), (1, 0), (24, 0), (1, 1))),
 ((1298, 1326), (31546951, 32210300, 0.663, 0), ((5, 1), (14, 1), (5, 0))),
 ((1326, 1371), (32210300, 32683165, 0.473, 0), ((0, 1), (9, 0), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (23, 0), (14, 0), (23, 1), (24, 1), (1, 0), (24, 0), (1, 1))),
 ((1326, 1371), (32210300, 32683165, 0.473, 0), ((5, 1), (14, 1), (5, 0))),
 ((1371, 1385), (32683165, 32850349, 0.167, 0), ((0, 1), (9, 0), (9, 1), (26, 1), (7, 1), (26, 0), (23, 0), (14, 0), (23, 1), (24, 1), (24, 0), (1, 1))),
 ((1371, 1385), (32683165, 32850349, 0.167, 0), ((5, 1), (14, 1), (5, 0))),
 ((1385, 1389), (32850349, 32878244, 0.028, 0), ((0, 1), (9, 1), (26, 1), (7, 1), (23, 0), (14, 0), (23, 1), (24, 1), (24, 0), (1, 1))),
 ((1385, 1389), (32850349, 32878244, 0.028, 0), ((5, 1), (14, 1), (5, 0))),
 ((1389, 1551), (32878244, 33995211, 1.117, 0), ((0, 1), (9, 1), (26, 1), (7, 1), (14, 0), (23, 1), (24, 1), (1, 1))),
 ((1389, 1551), (32878244, 33995211, 1.117, 0), ((14, 1), (5, 0))),
 ((1551, 1632), (33995211, 34523488, 0.528, 0), ((0, 1), (9, 1), (26, 1), (7, 1), (26, 0), (14, 0), (23, 1), (24, 1), (1, 1))),
 ((1551, 1632), (33995211, 34523488, 0.528, 0), ((14, 1), (5, 0))),
 ((1632, 1798), (34523488, 35412768, 0.889, 0), ((0, 1), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (14, 0), (23, 1), (24, 1), (1, 1))),
 ((1632, 1798), (34523488, 35412768, 0.889, 0), ((14, 1), (5, 0))),
 ((1798, 1802), (35412768, 35440705, 0.028, 0), ((0, 1), (7, 0), (9, 1), (26, 1), (7, 1), (26, 0), (23, 0), (14, 0), (23, 1), (24, 1), (24, 0), (1, 1))),
 ((1798, 1802), (35412768, 35440705, 0.028, 0), ((5, 1), (14, 1), (5, 0))),
 ((1802, 1864), (35440705, 35893156, 0.452, 0), ((0, 1), (7, 0), (9, 1), (26, 1), (7, 1), (23, 0), (14, 0), (23, 1), (24, 1), (24, 0), (1, 1))),
 ((1802, 1864), (35440705, 35893156, 0.452, 0), ((5, 1), (14, 1), (5, 0))),
 ((1864, 1948), (35893156, 36934821, 1.042, 0), ((0, 1), (23, 0), (14, 0), (26, 1), (23, 1), (24, 1), (24, 0))),
 ((1864, 1948), (35893156, 36934821, 1.042, 0), ((5, 1), (14, 1), (5, 0))),
 ((1948, 1986), (36934821, 37311806, 0.377, 0), ((0, 1), (0, 0), (26, 1), (23, 0), (14, 0), (23, 1), (24, 1), (24, 0))),
 ((1948, 1986), (36934821, 37311806, 0.377, 0), ((5, 1), (14, 1), (5, 0))),
 ((1986, 1994), (37311806, 37323988, 0.012, 0), ((0, 1), (0, 0), (14, 0), (26, 1), (23, 1), (24, 1))),
 ((1986, 1994), (37311806, 37323988, 0.012, 0), ((14, 1), (5, 0))),
 ((1994, 1998), (37323988, 37338286, 0.014, 0), ((0, 1), (0, 0), (14, 0), (26, 1), (23, 1), (24, 1), (26, 0))),
 ((1994, 1998), (37323988, 37338286, 0.014, 0), ((14, 1), (5, 0))),
 ((1998, 2000), (37338286, 37367308, 0.029, 0), ((0, 1), (0, 0), (26, 1), (26, 0), (23, 0), (14, 0), (23, 1), (24, 1), (24, 0))),
 ((1998, 2000), (37338286, 37367308, 0.029, 0), ((5, 1), (14, 1), (5, 0))),
 ((2000, 2004), (37367308, 37375117, 0.008, 0), ((0, 1), (0, 0), (26, 0), (23, 0), (14, 0), (23, 1), (24, 1), (24, 0))),
 ((2000, 2004), (37367308, 37375117, 0.008, 0), ((5, 1), (14, 1), (5, 0))),
 ((2004, 2087), (37375117, 38818676, 1.444, 0), ((0, 1), (0, 0), (26, 1), (26, 0), (23, 0), (14, 0), (23, 1), (24, 1), (24, 0))),
 ((2004, 2087), (37375117, 38818676, 1.444, 0), ((5, 1), (14, 1), (5, 0))),
 ((2087, 2093), (38818676, 38907849, 0.089, 0), ((26, 0), (14, 0), (26, 1))),
 ((2093, 2121), (38907849, 39197939, 0.290, 0), ((14, 0), (26, 1))),
 ((2129, 2219), (39320828, 41090009, 1.769, 0), ((26, 0), (14, 0))),
 ((2233, 2298), (41369108, 42654327, 1.285, 0), ((14, 0), (26, 1))),
 ((2298, 2462), (42654327, 44391234, 1.737, 0), ((26, 0), (14, 0), (26, 1))),
 ((2462, 2483), (44391234, 44557063, 0.166, 0), ((26, 0), (14, 0))),
 ((2483, 3108), (44557063, 49669972, 5.113, 0), ((26, 0), (14, 0), (26, 1))),
 ((3108, 3218), (49669972, 51156934, 1.487, 0), ((26, 0), (14, 0), (26, 1))),
 ((3108, 3218), (49669972, 51156934, 1.487, 0), ((14, 1), (6, 1)))],
                                     decimal=3, err_msg='Wrong IBD segments') 

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------