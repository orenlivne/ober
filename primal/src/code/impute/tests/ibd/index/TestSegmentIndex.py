'''
============================================================
Test IBD segment index queries.

Created on January 31, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, impute as im
#from numpy.testing.utils import assert_equal
#from numpy.ma.testutils import assert_almost_equal
#from nose.tools import nottest

class TestSegmentIndex(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------       
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.ibd = im.index.segment_index.SegmentIndex(im.itu.SEGMENT_INDEX_CHR22)
        self.snp = 1200
        
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
#    @nottest # Interface to SegmentIndex changed. Need to update index segment testdata files at some point         
#    def test_find_ibd_group_boundary_cast(self):
#        '''Test finding the IBD group of a haplotype in the last IBD group.
#        This caused index out-of-bounds error on 25-MAR-2013.'''
#        group = self.ibd.find(22, self.snp, 244, 0)
#        assert_equal(group, [[244, 1], [244, 0]], 'Wrong IBD group')
#        
#    @nottest  # Interface to SegmentIndex changed. Need to update index segment testdata files at some point         
#    def test_covered_by_training(self):
#        '''Test that % haps covered by WGS haps is correct. Predicts the allele call rate
#        if all WGS samples are phased.'''
#        covered = (1.0 * len(self.ibd.covered_by_training(22, self.snp, im.wgs_sample_index()))) / (2 * self.ibd.num_samples)
#        assert_almost_equal(covered, 0.9074, 4, 'Incorrect allele call rate')
#         
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
