'''
============================================================
Test IBD operations that do not fall under a specific IBD
implementation category (GERMLINE, IBDLD, etc.).

Created on July 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np
import impute.tools.genotype_tools as gt
from impute.ibd import ibd_family as ip
from impute import impute_test_util as itu
from impute.ibd import ibd, diff 
from numpy.ma.testutils import assert_equal
from impute.impute_test_util import Templates, assert_segments_almost_equal
from impute.tools.param import PhaseParam

class TestIbd(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------

    def setUp(self):
        '''Load single nuclear family test case.'''
        self.problem = Templates.problem_family(itu.FAMILY7, haplotype=True)
        self.haplotype = self.problem.haplotype
        self.family = self.problem.families(3)[0]
        (self.father, self.mother) = (self.family.father, self.family.mother)
        self.children = self.family.children
        
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_hap_diff(self):
        '''Check haplotype difference functionality.'''
        a = np.array([2, 2, 2, 1, 2, 1, 1, 0], dtype=np.byte)
        b = np.array([0, 2, 2, 2, 1, 1, 1, 2], dtype=np.byte)
        d = np.array([-1, 0, 0, 1, 1, 0, 0, -1], dtype=np.byte)
        assert_equal(diff.hap_diff(a, b), d, 'Wrong haplotype difference')
        
    def test_ibs_segments(self):
        '''Test locating IBD segments between parent and child haplotypes.'''
        (parent, child, hap1, hap2) = (1, 2, 0, 1)
        het_snps = gt.where_heterozygous(self.haplotype.data, parent)
        
        segment_set = ibd.ibs_segments(self.haplotype, parent, child, hap1, hap2, het_snps, True)
        assert_segments_almost_equal(segment_set,
                                     [((3   , 172), (17080378, 18541497, 1.461, 0), ((1, 1), (2, 1))),
                                      ((193 , 3218), (18996226, 51156934, 32.161, 1), ((1, 0), (2, 1)))],
                                     full_data=True, decimal=3, err_msg='Wrong IBD segments') 

        segment_set = ibd.ibs_segments(self.haplotype, parent, child, hap1, hap2, het_snps, True,
                                       length_bound='base_pair', min_segment_length=10.0)
        assert_segments_almost_equal(segment_set,
                                     [((193 , 3218), (18996226, 51156934, 32.161, 1), ((1, 0), (2, 1)))],
                                     full_data=True, decimal=3, err_msg='Wrong IBD segments') 

    def test_ibd_segments_in_family(self):
        '''Test locating recombination locations and subsequently IBD regions between 
        parent and child haplotypes.'''
        segment_set = ibd.concatenate_segments(ip.ibd_segments_in_family(self.haplotype, self.family,
                                                                         PhaseParam.HET_FILL_THRESHOLD))
        assert_segments_almost_equal(segment_set,
                                   [((3   , 172), (17080378, 18541497, 1.461, 0), ((1, 1), (2, 1))),
                                    ((193 , 3218), (18996226, 51156934, 32.161, 1), ((1, 0), (2, 1))),
                                    ((0   , 3218), (16484792, 51156934, 34.672, 4), ((1, 0), (3, 1))),
                                    ((3   , 240), (17080378, 19593301, 2.513, 0), ((1, 0), (4, 1))),
                                    ((256 , 3218), (19873357, 51156934, 31.284, 0), ((4, 1), (1, 1))),
                                    ((0   , 3218), (16484792, 51156934, 34.672, 0), ((5, 1), (1, 0))),
                                    ((3   , 2696), (17080378, 46630634, 29.550, 0), ((6, 1), (1, 1))),
                                    ((2718, 3218), (46853292, 51156934, 4.304, 1), ((1, 0), (6, 1))),
                                    ((0   , 3218), (16484792, 51156934, 34.672, 0), ((1, 0), (7, 1))),
                                    ((3   , 1908), (17080378, 36572515, 19.492, 0), ((8, 1), (1, 0))),
                                    ((1939, 3218), (36804282, 51156934, 14.353, 0), ((8, 1), (1, 1)))],
                                    full_data=True, decimal=3, err_msg='Wrong IBD segments')
        assert_equal(segment_set.errors,
                     [[2997, 2997, 132, 1364, 2493, 2800, 132, 1364, 2493, 2800, 2997, 2997],
                      [   1, 2, 1, 1, 1, 1, 3, 3, 3, 3, 1, 6]],
                     'Wrong genotype errors')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
