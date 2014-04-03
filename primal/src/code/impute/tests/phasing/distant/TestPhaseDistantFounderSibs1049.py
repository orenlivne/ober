'''
============================================================
Test distant-IBD (segments IBS-based for now) phasing for
unphased individuals in a family of founder sibs, i.e.,
genotyped sibs whose parents are not genotyped.  

In the considered sib founder family (family of 507),
the unphased individual is 1049. In the sub-problem indices,
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
from impute.tools.param import PhaseParam

class TestPhaseDistantFounderSibs1049(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        self.problem = im.io.read_npz(im.itu.SIB_FOUNDERS_1049_STAGE3)
        self.phaser = im.phase_core.new_phaser_chain([im.phase_distant.family_sib_comparison_phaser()])
        self.unphased_samples = self.problem.find_samples_with_fill_lt(0.9, sample=self.problem.first_family.member_list)
        self.s = int(self.unphased_samples[0][0])
        self.params = im.param.PhaseParam()
        self.proband = 3
        
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
        assert_almost_equal(self.unphased_samples, [[ self.proband, 0.63424487]],
                            decimal=5, err_msg='Wrong unphased samples')

    def test_find_phased_sibs(self):
        '''Test locating phased sibs of the unphased sample.'''
        relatives = self.__phased_sibs(self.s)
        assert_equal(relatives, [0, 1, 2, 4], err_msg='Wrong phased sib set - IDs')
        for sample in relatives:
            self.assertTrue(self.problem.haplotype.fill_fraction(sample=sample) >= self.params.surrogate_parent_fill_threshold,
                            'Surrogate parent''s not phased enough')
        
    def test_family_sib_comparison_phaser(self):
        '''Check phasing the test data set up to stage 3 (nuclear families). The fill % will
        be low here since we're dealing with non-nuclear-family-members.'''
        im.itu.assert_size_equals(self.problem.genotype, 3218, 5)
        im.itu.assert_problem_stats(self.problem, 32180, 29597, 70)
        h = self.problem.haplotype
        assert_almost_equal(h.fill_fraction(sample=self.proband), 0.634, 3, 'Unexpected pre-phasing fill %')
        self.phaser.run(self.problem)
        im.itu.assert_problem_stats(self.problem, 32180, 31946, 70)
        assert_almost_equal(h.fill_fraction(sample=self.proband), 0.993, 3, 'Unexpected post-phasing fill %')

    def test_ibd_segments_hmm(self):
        '''Test locating IBD segments between the unphased proband and its sib surrogate parents.
        Uses HMM IBD posterior.'''
        relatives = self.__phased_sibs(self.s)
        segment_set = im.idist.ibd_segments_with_relatives(self.problem, self.s, relatives, PhaseParam(debug=False),
                                                           im.ibd_hmm.prob_ibd_hmm)
        segment_set.group_to_disjoint()
        assert_segments_almost_equal(segment_set,
[((0   , 96  ), (16484792, 17948473,   1.464, 0), ((2, 0),(3, 1),(4, 1),(2, 1))),
 ((96  , 344 ), (17948473, 21460008,   3.512, 0), ((3, 1),(4, 1),(2, 1))),
 ((344 , 380 ), (21460008, 22554306,   1.094, 0), ((2, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((380 , 383 ), (22554306, 22555078,   0.001, 0), ((3, 1),(4, 1),(2, 1))),
 ((383 , 438 ), (22555078, 23636541,   1.081, 0), ((2, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((438 , 442 ), (23636541, 23695404,   0.059, 0), ((3, 1),(4, 1),(2, 1))),
 ((442 , 519 ), (23695404, 24406778,   0.711, 0), ((2, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((519 , 557 ), (24406778, 25088629,   0.682, 0), ((0, 1),(2, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((557 , 951 ), (25088629, 27698217,   2.610, 0), ((0, 1),(3, 1),(4, 1),(2, 1))),
 ((951 , 969 ), (27698217, 27832985,   0.135, 0), ((0, 1),(2, 0),(3, 1),(4, 1),(2, 1))),
 ((969 , 1019), (27832985, 28093392,   0.260, 0), ((2, 0),(3, 1),(4, 1),(2, 1))),
 ((1019, 1147), (28093392, 29670939,   1.578, 0), ((2, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((1147, 1188), (29670939, 30113960,   0.443, 0), ((2, 0),(3, 1),(4, 1),(2, 1))),
 ((1188, 1403), (30113960, 32950053,   2.836, 0), ((2, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((1403, 1943), (32950053, 36891858,   3.942, 0), ((2, 0),(3, 1),(4, 1),(2, 1))),
 ((1943, 2053), (36891858, 37982012,   1.090, 0), ((2, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((2053, 2055), (37982012, 38086574,   0.105, 0), ((2, 0),(3, 1),(4, 1),(2, 1))),
 ((2055, 2133), (38086574, 39454432,   1.368, 0), ((2, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((2133, 2174), (39454432, 40018212,   0.564, 0), ((2, 0),(3, 1),(4, 1),(2, 1))),
 ((2174, 2221), (40018212, 41107688,   1.089, 0), ((2, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((2221, 2612), (41107688, 45515269,   4.408, 0), ((2, 0),(3, 1),(4, 1),(2, 1))),
 ((2612, 2661), (45515269, 45972017,   0.457, 0), ((2, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((2661, 2735), (45972017, 47094390,   1.122, 0), ((0, 1),(0, 0),(3, 1),(2, 1),(2, 0),(1, 0),(4, 1),(1, 1),(4, 0))),
 ((2735, 2945), (47094390, 48569604,   1.475, 0), ((2, 0),(0, 0),(1, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((2945, 2991), (48569604, 48741583,   0.172, 0), ((2, 0),(0, 0),(1, 0),(3, 1),(4, 1),(4, 0))),
 ((2991, 3127), (48741583, 49752332,   1.011, 0), ((2, 0),(0, 0),(1, 0),(3, 1),(4, 1),(2, 1),(4, 0))),
 ((3127, 3177), (49752332, 50120255,   0.368, 0), ((2, 0),(1, 0),(3, 1),(4, 1),(0, 0))),
 ((3177, 3218), (50120255, 51156934,   1.037, 0), ((0, 1),(0, 0),(3, 1),(2, 1),(2, 0),(1, 0),(4, 1),(1, 1),(4, 0)))],
                                      decimal=3, err_msg='Wrong IBD segments') 

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
