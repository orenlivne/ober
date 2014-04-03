'''
============================================================
Test phasing within a nuclear family that Gaixin found to
possibly have too much zeroed-out genotypes because of
disagreement among the children that could be due to an
error in the parent only.

Created on October 19, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest
# import numpy as np
from impute.ibd import ibd_family as ip
from impute import impute_test_util as itu
from numpy.ma.testutils import assert_equal
from impute.phasing.phase_trivial import trivial_phaser
from impute.tools import genotype_tools as gt
from impute.data import io
from impute.phasing import phase_core
from impute.ibd import ibd, ibd_child as ic 
from chain import FilterChain
from impute.phasing.phase_family import family_phaser, family_child_comparison_phaser
from impute.data.constants import MATERNAL
from impute.impute_test_util import assert_segments_almost_equal
from impute.tools.param import PhaseParam
from util import Struct
import itertools
from impute.ibd.segment import SegmentSet, Segment

class TestPhaseFamilyCompareChild0(unittest.TestCase):
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
        self.problem = io.read_npz(itu.FAMILY_TOO_ZEROED_STAGE1)
        self.family = self.problem.families()[0]
        self.phaser = phase_core.PhaseDecorator(FilterChain([trivial_phaser(),
                                                             family_phaser(),
                                                             family_child_comparison_phaser()]))
        self.comparator = ic.ChildComparator(Struct(problem=self.problem, params=PhaseParam()), self.family)

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_ibd_segments_parent_child(self):
        '''Test IBD segments that are apparently cut too short between a parent and child.'''
        parent = self.family.mother
        child = 2
        parent_type = MATERNAL
        h = self.problem.haplotype
        het_snps = gt.where_heterozygous(h.data, parent)
        segment_set = ip.ibd_segments_parent_child(h, parent, child, parent_type, het_snps)
        assert_segments_almost_equal(segment_set,
                                     [((2   , 2650), (17075353, 45892433, 28.817, 0), ((1, 1), (2, 1))),
                                      ((2657, 3218), (45940934, 51156934, 5.216, 2), ((1, 0), (2, 1)))],
                                    full_data=True, decimal=3, err_msg='Wrong IBD segments')

    def test_ibd_segments_in_family(self):
        '''Test locating recombination locations and subsequently IBD regions between 
        parent and child haplotypes.'''
        segment_set = ibd.concatenate_segments(ip.ibd_segments_in_family(self.problem.haplotype, self.family,
                                                                         PhaseParam.HET_FILL_THRESHOLD))
        assert_segments_almost_equal(segment_set,
                                   [((4   , 1411), (17087656, 32978753, 15.891, 0), ((2, 0), (0, 0))),
                                    ((1413, 3215), (33013062, 51103692, 18.091, 1), ((0, 1), (2, 0))),
                                    ((0   , 3218), (16484792, 51156934, 34.672, 1), ((0, 1), (3, 0))),
                                    ((4   , 2551), (17087656, 45007348, 27.920, 1), ((0, 1), (4, 0))),
                                    ((2569, 3215), (45198494, 51103692, 5.905, 0), ((0, 0), (4, 0))),
                                    ((0   , 3218), (16484792, 51156934, 34.672, 1), ((0, 0), (5, 0))),
                                    ((4   , 795), (17087656, 27011420, 9.924, 0), ((0, 0), (6, 0))),
                                    ((805 , 3215), (27119061, 51103692, 23.985, 1), ((0, 1), (6, 0))),
                                    ((2   , 2650), (17075353, 45892433, 28.817, 0), ((1, 1), (2, 1))),
                                    ((2657, 3218), (45940934, 51156934, 5.216, 2), ((1, 0), (2, 1))),
                                    ((1   , 576), (17065079, 25228375, 8.163, 0), ((3, 1), (1, 1))),
                                    ((600 , 3218), (25444874, 51156934, 25.712, 2), ((1, 0), (3, 1))),
                                    ((2   , 1978), (17075353, 37228277, 20.153, 2), ((1, 0), (4, 1))),
                                    ((2019, 3218), (37509844, 51156934, 13.647, 0), ((4, 1), (1, 1))),
                                    ((0   , 3218), (16484792, 51156934, 34.672, 3), ((5, 1), (1, 0))),
                                    ((2   , 301), (17075353, 20993766, 3.918, 0), ((6, 1), (1, 1))),
                                    ((334 , 3218), (21363960, 51156934, 29.793, 2), ((1, 0), (6, 1)))],
                                    full_data=True, decimal=3, err_msg='Wrong IBD segments')

        # Check that there are at most 4 distinct haplotypes at each SNP        
        segment_set.group_to_disjoint(True)
        assert_segments_almost_equal(segment_set,
                                     [((0   , 1), (16484792, 17065079, 0.580, 0), ((5, 1), (1, 0))),
                                      ((0   , 1), (16484792, 17065079, 0.580, 0), ((0, 1), (3, 0))),
                                      ((0   , 1), (16484792, 17065079, 0.580, 0), ((0, 0), (5, 0))),
                                      ((1   , 2), (17065079, 17075353, 0.010, 0), ((5, 1), (1, 0))),
                                      ((1   , 2), (17065079, 17075353, 0.010, 0), ((3, 1), (1, 1))),
                                      ((1   , 2), (17065079, 17075353, 0.010, 0), ((0, 1), (3, 0))),
                                      ((1   , 2), (17065079, 17075353, 0.010, 0), ((0, 0), (5, 0))),
                                      ((2   , 4), (17075353, 17087656, 0.012, 0), ((0, 0), (5, 0))),
                                      ((2   , 4), (17075353, 17087656, 0.012, 0), ((0, 1), (3, 0))),
                                      ((2   , 4), (17075353, 17087656, 0.012, 0), ((6, 1), (3, 1), (1, 1), (2, 1))),
                                      ((2   , 4), (17075353, 17087656, 0.012, 0), ((5, 1), (1, 0), (4, 1))),
                                      ((4   , 301), (17087656, 20993766, 3.906, 0), ((2, 0), (0, 0), (6, 0), (5, 0))),
                                      ((4   , 301), (17087656, 20993766, 3.906, 0), ((0, 1), (3, 0), (4, 0))),
                                      ((4   , 301), (17087656, 20993766, 3.906, 0), ((6, 1), (3, 1), (1, 1), (2, 1))),
                                      ((4   , 301), (17087656, 20993766, 3.906, 0), ((5, 1), (1, 0), (4, 1))),
                                      ((301 , 334), (20993766, 21363960, 0.370, 0), ((2, 0), (0, 0), (6, 0), (5, 0))),
                                      ((301 , 334), (20993766, 21363960, 0.370, 0), ((0, 1), (3, 0), (4, 0))),
                                      ((301 , 334), (20993766, 21363960, 0.370, 0), ((3, 1), (1, 1), (2, 1))),
                                      ((301 , 334), (20993766, 21363960, 0.370, 0), ((5, 1), (1, 0), (4, 1))),
                                      ((334 , 576), (21363960, 25228375, 3.864, 0), ((2, 0), (0, 0), (6, 0), (5, 0))),
                                      ((334 , 576), (21363960, 25228375, 3.864, 0), ((0, 1), (3, 0), (4, 0))),
                                      ((334 , 576), (21363960, 25228375, 3.864, 0), ((3, 1), (1, 1), (2, 1))),
                                      ((334 , 576), (21363960, 25228375, 3.864, 0), ((5, 1), (6, 1), (4, 1), (1, 0))),
                                      ((576 , 600), (25228375, 25444874, 0.216, 0), ((1, 1), (2, 1))),
                                      ((576 , 600), (25228375, 25444874, 0.216, 0), ((2, 0), (0, 0), (6, 0), (5, 0))),
                                     ((576 , 600), (25228375, 25444874, 0.216, 0), ((0, 1), (3, 0), (4, 0))),
                                     ((576 , 600), (25228375, 25444874, 0.216, 0), ((5, 1), (6, 1), (4, 1), (1, 0))),
                                     ((600 , 795), (25444874, 27011420, 1.567, 0), ((5, 1), (6, 1), (3, 1), (4, 1), (1, 0))),
                                     ((600 , 795), (25444874, 27011420, 1.567, 0), ((2, 0), (0, 0), (6, 0), (5, 0))),
                                     ((600 , 795), (25444874, 27011420, 1.567, 0), ((0, 1), (3, 0), (4, 0))),
                                     ((600 , 795), (25444874, 27011420, 1.567, 0), ((1, 1), (2, 1))),
                                     ((795 , 805), (27011420, 27119061, 0.108, 0), ((5, 1), (6, 1), (3, 1), (4, 1), (1, 0))),
                                     ((795 , 805), (27011420, 27119061, 0.108, 0), ((2, 0), (0, 0), (5, 0))),
                                     ((795 , 805), (27011420, 27119061, 0.108, 0), ((0, 1), (3, 0), (4, 0))),
                                     ((795 , 805), (27011420, 27119061, 0.108, 0), ((1, 1), (2, 1))),
                                     ((805 , 1411), (27119061, 32978753, 5.860, 0), ((5, 1), (6, 1), (3, 1), (4, 1), (1, 0))),
                                     ((805 , 1411), (27119061, 32978753, 5.860, 0), ((2, 0), (0, 0), (5, 0))),
                                     ((805 , 1411), (27119061, 32978753, 5.860, 0), ((0, 1), (3, 0), (6, 0), (4, 0))),
                                     ((805 , 1411), (27119061, 32978753, 5.860, 0), ((1, 1), (2, 1))),
                                     ((1411, 1413), (32978753, 33013062, 0.034, 0), ((5, 1), (6, 1), (3, 1), (4, 1), (1, 0))),
                                     ((1411, 1413), (32978753, 33013062, 0.034, 0), ((1, 1), (2, 1))),
                                     ((1411, 1413), (32978753, 33013062, 0.034, 0), ((0, 1), (3, 0), (6, 0), (4, 0))),
                                     ((1411, 1413), (32978753, 33013062, 0.034, 0), ((0, 0), (5, 0))),
                                     ((1413, 1978), (33013062, 37228277, 4.215, 0), ((5, 1), (6, 1), (3, 1), (4, 1), (1, 0))),
                                     ((1413, 1978), (33013062, 37228277, 4.215, 0), ((0, 1), (3, 0), (4, 0), (6, 0), (2, 0))),
                                     ((1413, 1978), (33013062, 37228277, 4.215, 0), ((1, 1), (2, 1))),
                                     ((1413, 1978), (33013062, 37228277, 4.215, 0), ((0, 0), (5, 0))),
                                     ((1978, 2019), (37228277, 37509844, 0.282, 0), ((1, 1), (2, 1))),
                                     ((1978, 2019), (37228277, 37509844, 0.282, 0), ((0, 1), (3, 0), (4, 0), (6, 0), (2, 0))),
                                     ((1978, 2019), (37228277, 37509844, 0.282, 0), ((5, 1), (6, 1), (3, 1), (1, 0))),
                                     ((1978, 2019), (37228277, 37509844, 0.282, 0), ((0, 0), (5, 0))),
                                     ((2019, 2551), (37509844, 45007348, 7.498, 0), ((0, 1), (3, 0), (4, 0), (6, 0), (2, 0))),
                                     ((2019, 2551), (37509844, 45007348, 7.498, 0), ((0, 0), (5, 0))),
                                     ((2019, 2551), (37509844, 45007348, 7.498, 0), ((5, 1), (6, 1), (3, 1), (1, 0))),
                                     ((2019, 2551), (37509844, 45007348, 7.498, 0), ((4, 1), (1, 1), (2, 1))),
                                     ((2551, 2569), (45007348, 45198494, 0.191, 0), ((0, 0), (5, 0))),
                                     ((2551, 2569), (45007348, 45198494, 0.191, 0), ((0, 1), (3, 0), (6, 0), (2, 0))),
                                     ((2551, 2569), (45007348, 45198494, 0.191, 0), ((5, 1), (6, 1), (3, 1), (1, 0))),
                                     ((2551, 2569), (45007348, 45198494, 0.191, 0), ((4, 1), (1, 1), (2, 1))),
                                     ((2569, 2650), (45198494, 45892433, 0.694, 0), ((0, 0), (5, 0), (4, 0))),
                                     ((2569, 2650), (45198494, 45892433, 0.694, 0), ((0, 1), (3, 0), (6, 0), (2, 0))),
                                     ((2569, 2650), (45198494, 45892433, 0.694, 0), ((5, 1), (6, 1), (3, 1), (1, 0))),
                                     ((2569, 2650), (45198494, 45892433, 0.694, 0), ((4, 1), (1, 1), (2, 1))),
                                     ((2650, 2657), (45892433, 45940934, 0.049, 0), ((0, 0), (5, 0), (4, 0))),
                                     ((2650, 2657), (45892433, 45940934, 0.049, 0), ((0, 1), (3, 0), (6, 0), (2, 0))),
                                     ((2650, 2657), (45892433, 45940934, 0.049, 0), ((5, 1), (6, 1), (3, 1), (1, 0))),
                                     ((2650, 2657), (45892433, 45940934, 0.049, 0), ((4, 1), (1, 1))),
                                     ((2657, 3215), (45940934, 51103692, 5.163, 0), ((5, 1), (6, 1), (3, 1), (1, 0), (2, 1))),
                                     ((2657, 3215), (45940934, 51103692, 5.163, 0), ((0, 1), (3, 0), (6, 0), (2, 0))),
                                     ((2657, 3215), (45940934, 51103692, 5.163, 0), ((0, 0), (5, 0), (4, 0))),
                                     ((2657, 3215), (45940934, 51103692, 5.163, 0), ((4, 1), (1, 1))),
                                     ((3215, 3218), (51103692, 51156934, 0.053, 0), ((5, 1), (6, 1), (3, 1), (1, 0), (2, 1))),
                                     ((3215, 3218), (51103692, 51156934, 0.053, 0), ((4, 1), (1, 1))),
                                     ((3215, 3218), (51103692, 51156934, 0.053, 0), ((0, 1), (3, 0))),
                                     ((3215, 3218), (51103692, 51156934, 0.053, 0), ((0, 0), (5, 0)))],
                                     full_data=True, decimal=3, err_msg='Wrong IBD segments')
        assert_equal(segment_set.errors,
                     [[2162, 2162, 2162, 2162, 2162, 2162, 2162, 2162, 2162, 2162, 2846, 3202, 2846,
                       3202, 1092, 3202, 1092, 3202, 3, 1092, 3, 1092, 3, 1092, 3202, 3, 1092, 3202,
                        1092, 3202, 1092, 3202],
                      [0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 1, 1, 2, 2, 1, 1, 3, 3, 1, 1, 4, 4, 1, 1, 1, 5, 5,
                       5, 1, 1, 6, 6]],
                     'Wrong genotype errors')
        
        self.assertTrue(max(len(list(g)) for (_, g) in itertools.groupby(segment_set, lambda segment: segment.snp)) <= 4,
                        'Too many distinct haplotypes in nuclear family')
        
        # Extract IBD segments of child's haplotype vs. sibs
        segments_of_30 = SegmentSet(Segment(x.snp, x.samples - set([(0, 0), (0, 1), (1, 0), (1, 1)]), x.bp, error_snps=x.error_snps)
                                   for x in segment_set if (3, 0) in x.samples)
        segments_of_30.merge_consecutive()
        assert_segments_almost_equal(segments_of_30,
                                     [((0   , 4), (16484792, 17087656, 0.603, 0), ((3, 0),)),
                                      ((4   , 805), (17087656, 27119061, 10.031, 0), ((3, 0), (4, 0))),
                                      ((805 , 1413), (27119061, 33013062, 5.894, 0), ((3, 0), (6, 0), (4, 0))),
                                      ((1413, 2551), (33013062, 45007348, 11.994, 0), ((3, 0), (2, 0), (6, 0), (4, 0))),
                                      ((2551, 3215), (45007348, 51103692, 6.096, 0), ((3, 0), (2, 0), (6, 0))),
                                      ((3215, 3218), (51103692, 51156934, 0.053, 0), ((3, 0),))],
                                     full_data=True, decimal=3, err_msg='Wrong IBD segments')

        segments_of_31 = SegmentSet(Segment(x.snp, x.samples - set([(0, 0), (0, 1), (1, 0), (1, 1)]), x.bp, error_snps=x.error_snps)
                                   for x in segment_set if (3, 1) in x.samples)
        segments_of_31.merge_consecutive()
        assert_segments_almost_equal(segments_of_31,
                                     [((1   , 2), (17065079, 17075353, 0.010, 0), ((3, 1),)),
                                      ((2   , 301), (17075353, 20993766, 3.918, 0), ((6, 1), (3, 1), (2, 1))),
                                      ((301 , 576), (20993766, 25228375, 4.235, 0), ((3, 1), (2, 1))),
                                      ((600 , 1978), (25444874, 37228277, 11.783, 0), ((5, 1), (6, 1), (3, 1), (4, 1))),
                                      ((1978, 2657), (37228277, 45940934, 8.713, 0), ((5, 1), (6, 1), (3, 1))),
                                      ((2657, 3218), (45940934, 51156934, 5.216, 0), ((5, 1), (6, 1), (3, 1), (2, 1)))],
                                     full_data=True, decimal=3, err_msg='Wrong IBD segments')
             
    def test_outer_duo(self):
        '''Test applying child comparison to a nuclear family with many genotyped kids but only
        one genotyped parent. Seems to be fine for now: too many errors are flagged, but we are not
        going to split hair.'''
        p = self.problem
        # h = p.haplotype
        # (f, m) = (self.family.father, self.family.mother)
        snp = p.info.snp_by_name('rs5746679')
        assert_equal(snp, [3], 'Wrong SNP index')
        itu.assert_size_equals(p.genotype, 3218, 7)

        itu.assert_problem_stats(p, 45052, 40086, 4)
#        print 'genotypes'
#        print p.genotype.data[snp, :, :]
#        print 'haplotypes'
#        print p.haplotype.data[snp, :, :]

        phaser = family_phaser()
        phaser.run(p)
        
        itu.assert_problem_stats(p, 45052, 45023, 25)

#        print 'genotypes'
#        print p.genotype.data[snp, :, :]
#        print 'haplotypes'
#        print p.haplotype.data[snp, :, :]

#    def test_child_comparison(self):
#        '''Test applying child comparison to a nuclear family with many genotyped kids but only
#        one genotyped parent.'''
#        p = self.problem
#        h = p.haplotype
#        (f, m) = (self.family.father, self.family.mother)
#        itu.assert_size_equals(p.genotype, 3218, 7)
#        itu.assert_problem_stats(p, 45052, 45023, 25)
#        phaser = family_child_comparison_phaser(debug=False)
#        phaser.run(p)
#        itu.assert_problem_stats(p, 45052, 45023, 25)
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
