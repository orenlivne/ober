'''
============================================================
Test distant-IBD (segments IBS-based for now) phasing for
unphased individuals in a family of founder sibs, i.e.,
genotyped sibs whose parents are not genotyped.  

In the considered sib founder family (family of 507),
the unphased individual is 139. In the sub-problem indices,
it is sample #2. Fill fraction before distant IBD:

p.fill_fraction() =
array([[ 0.        ,  0.99720323],
       [ 1.        ,  0.99829086],
       [ 2.        ,  0.63704164], <-----
       [ 3.        ,  0.99627098],
       [ 4.        ,  0.99704786],
       [ 5.        ,  0.9883468 ]])


Created on December 10, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, impute as im
from numpy.ma.testutils import assert_almost_equal, assert_equal
from impute.impute_test_util import assert_segments_almost_equal
#from impute.tools.param import PhaseParam

class TestPhaseDistantFounderSibs(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        self.problem = im.io.read_npz(im.itu.SIB_FOUNDERS_STAGE3)
        self.phaser = im.phase_core.new_phaser_chain([im.phase_distant.family_sib_comparison_phaser()])
        self.unphased_samples = self.problem.find_samples_with_fill_lt(0.9, sample=self.problem.first_family.member_list)
        self.s = int(self.unphased_samples[0][0])
        self.params = im.param.PhaseParam()
        
    def __phased_sibs(self, sample):
        genotyped_children = im.gt.genotyped_children(self.problem, self.problem.first_family)
        # Phased children = relatives array; at distance two from target sample
        # Note: these are all full sibs of the proband, not half sibs
        phased_children = [x for x in genotyped_children if self.problem.haplotype.fill_fraction(sample=x) >= self.params.surrogate_parent_fill_threshold]
        return phased_children
        
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_find_unphased_sample(self):
        '''Test locating unphased samples.'''
        assert_almost_equal(self.unphased_samples, [[ 2.        , 0.63704164]],
                            decimal=5, err_msg='Wrong unphased samples')

    def test_find_phased_sibs(self):
        '''Test locating phased sibs of the unphased sample.'''
        relatives = self.__phased_sibs(self.s)
        assert_equal(relatives, [0, 1, 3, 4, 5], err_msg='Wrong phased sib set - IDs')
        for sample in relatives:
            self.assertTrue(self.problem.haplotype.fill_fraction(sample=sample) >= self.params.surrogate_parent_fill_threshold,
                            'Surrogate parent''s not phased enough')
        
    def test_family_sib_comparison_phaser(self):
        '''Check phasing the test data set up to stage 3 (nuclear families). The fill % will
        be low here since we're dealing with non-nuclear-family-members.'''
        im.itu.assert_size_equals(self.problem.genotype, 3218, 6)
        im.itu.assert_problem_stats(self.problem, 38616, 36133, 65)
        h = self.problem.haplotype
        assert_almost_equal(h.fill_fraction(sample=2), 0.64, 2, 'Unexpected pre-phasing fill %')
        self.phaser.run(self.problem)
        im.itu.assert_problem_stats(self.problem, 38616, 38444, 65)
        assert_almost_equal(h.fill_fraction(sample=2), 0.996, 2, 'Unexpected post-phasing fill %')

    def test_ibd_segments_ibs(self):
        '''Test locating IBD segments between distant individuals. No IBD posterior filter (IBS=>IBD).'''
        relatives = self.__phased_sibs(self.s)
        segment_set = im.idist.ibd_segments_with_relatives(self.problem, self.s, relatives, self.params,
                                                           im.ibd_distant.prob_ibd_ibs)
        segment_set.group_to_disjoint()
        assert_segments_almost_equal(segment_set,
[((0   , 38  ), (16484792, 17561913,   1.077, 0), ((5, 1),(1, 0),(1, 1),(2, 1))),
 ((38  , 180 ), (17561913, 18562888,   1.001, 0), ((1, 0),(1, 1),(2, 1))),
 ((180 , 195 ), (18562888, 19030336,   0.467, 0), ((1, 0),(1, 1),(5, 0),(2, 1))),
 ((195 , 329 ), (19030336, 21242544,   2.212, 0), ((1, 0),(5, 0),(2, 1))),
 ((329 , 369 ), (21242544, 22485301,   1.243, 0), ((5, 1),(1, 0),(5, 0),(2, 1))),
 ((369 , 435 ), (22485301, 23627369,   1.142, 0), ((1, 0),(5, 0),(2, 1))),
 ((435 , 506 ), (23627369, 24235198,   0.608, 0), ((1, 0),(1, 1),(5, 0),(2, 1))),
 ((506 , 581 ), (24235198, 25304473,   1.069, 0), ((5, 1),(1, 0),(1, 1),(5, 0),(2, 1))),
 ((581 , 865 ), (25304473, 27370273,   2.066, 0), ((1, 0),(1, 1),(5, 0),(2, 1))),
 ((865 , 1032), (27370273, 28120707,   0.750, 0), ((1, 0),(2, 1),(1, 1),(5, 0),(4, 0))),
 ((1032, 1037), (28120707, 28165334,   0.045, 0), ((5, 0),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((1037, 1109), (28165334, 29400515,   1.235, 0), ((5, 0),(5, 1),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((1109, 1138), (29400515, 29495373,   0.095, 0), ((1, 0),(2, 1),(1, 1),(5, 0),(4, 0))),
 ((1138, 1249), (29495373, 30864692,   1.369, 0), ((5, 0),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((1249, 1276), (30864692, 31207207,   0.343, 0), ((1, 0),(2, 1),(1, 1),(5, 0),(4, 0))),
 ((1276, 1460), (31207207, 33291742,   2.085, 0), ((5, 0),(5, 1),(1, 0),(1, 1),(2, 1),(4, 0))),
 ((1460, 1792), (33291742, 35344317,   2.053, 0), ((1, 0),(2, 1),(1, 1),(5, 0),(4, 0))),
 ((1792, 1934), (35344317, 36746079,   1.402, 0), ((5, 0),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((1934, 2008), (36746079, 37383809,   0.638, 0), ((1, 0),(2, 1),(1, 1),(5, 0),(4, 0))),
 ((2008, 2052), (37383809, 37977282,   0.593, 0), ((5, 0),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((2052, 2068), (37977282, 38538246,   0.561, 0), ((5, 0),(5, 1),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((2068, 2170), (38538246, 39988277,   1.450, 0), ((5, 0),(5, 1),(1, 0),(1, 1),(2, 1),(4, 0))),
 ((2170, 2187), (39988277, 40351438,   0.363, 0), ((1, 0),(2, 1),(1, 1),(5, 0),(4, 0))),
 ((2187, 2189), (40351438, 40395084,   0.044, 0), ((5, 0),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((2189, 2269), (40395084, 42235352,   1.840, 0), ((5, 0),(5, 1),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((2269, 2313), (42235352, 42919125,   0.684, 0), ((5, 0),(5, 1),(1, 0),(1, 1),(2, 1),(4, 0))),
 ((2313, 2314), (42919125, 42982572,   0.063, 0), ((5, 0),(5, 1),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((2314, 2315), (42982572, 43013836,   0.031, 0), ((5, 0),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((2315, 2393), (43013836, 44023173,   1.009, 0), ((5, 0),(5, 1),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((2393, 2414), (44023173, 44101411,   0.078, 0), ((5, 0),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((2414, 2598), (44101411, 45409136,   1.308, 0), ((1, 0),(2, 1),(1, 1),(5, 0),(4, 0))),
 ((2598, 2661), (45409136, 45972017,   0.563, 0), ((1, 0),(1, 1),(2, 1),(4, 0))),
 ((2661, 2666), (45972017, 46060541,   0.089, 0), ((5, 1),(1, 0),(1, 1),(2, 1),(4, 0))),
 ((2666, 2744), (46060541, 47175815,   1.115, 0), ((5, 0),(5, 1),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((2744, 2779), (47175815, 47387316,   0.212, 0), ((5, 0),(5, 1),(1, 0),(1, 1),(2, 1),(4, 0))),
 ((2779, 2835), (47387316, 47786315,   0.399, 0), ((1, 0),(2, 1),(1, 1),(5, 0),(4, 0))),
 ((2835, 3043), (47786315, 49171641,   1.385, 0), ((1, 0),(1, 1),(2, 1),(4, 0))),
 ((3043, 3178), (49171641, 50120525,   0.949, 0), ((5, 0),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0))),
 ((3178, 3218), (50120525, 51156934,   1.036, 0), ((5, 0),(5, 1),(1, 0),(4, 1),(1, 1),(2, 1),(4, 0)))],
                                      decimal=3, err_msg='Wrong IBD segments') 

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
