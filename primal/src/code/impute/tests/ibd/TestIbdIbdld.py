'''
============================================================
Test the IBDLD interface for locating IBD segments.

Created on October 19, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, impute as im
#import impute.tools.genotype_tools as gt
#from impute.ibd import ibd_family as ip
from impute import impute_test_util as itu
from impute.ibd import segment#, ibd, diff 
from numpy.ma.testutils import assert_equal
from impute.impute_test_util import Templates, assert_segments_almost_equal
from impute.tools.genotype_tools import empty_errors_array
from impute.tools.param import PhaseParam

class TestIbdIbdld(unittest.TestCase):
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
    def test_ibd_segment_computer_ibdld(self):
        '''Test creating and reading from an IBDLD results cache. Segment unit is base pairs.'''
        cache = im.ibdld.ibd_ld.IbdSegmentGlobalCacheIbdld(itu.FAMILY7 + '.ibd')
        segments = cache.segments(22, 159351, 160541)
        assert_equal(segments, [[15445079, 17973301], [38435727, 41860926]],
                     'Wrong IBD segments read from cache')

    def test_ibd_segments_ibdld(self):
        '''Calculate IBD segments in a nuclear family using IBDLD.'''
        segment_cache = im.ibdld.ibd_ld.IbdSegmentGlobalCacheIbdld(itu.FAMILY7 + '.ibd')
        segment_computer = im.ibdld.ibd_ld.IbdSegmentComputerIbdld(segment_cache, self.haplotype,
                                                          chrom=22,
                                                          sample_id=self.problem.pedigree.sample_id,
                                                          samples=[2, 8],
                                                          threshold=0.9,
                                                          params=PhaseParam())
        segment_set = segment.break_and_group_segments(segment_computer.segments)
        assert_segments_almost_equal(segment_set,
                                     [((38  , 2849), [], ((8, 1), (2, 1)))],
                                     full_data=False, decimal=3, err_msg='Wrong grouped IBDLD IBD segments')
        assert_equal(segment_set.errors, empty_errors_array(), 'IBDLD does not support errors but they are output?!')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
