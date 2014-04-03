'''
============================================================
Test imputation routines.

Created on November 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, impute as im#, numpy as np
from numpy.ma.testutils import assert_equal#, assert_almost_equal
from impute.impute_test_util import assert_size_equals
from impute.imputation.ImputationSet import ImputationSet

class TestImputation(unittest.TestCase):
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

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_segments_intersecting_bp(self):
        '''Test finding all IBD dictionary segments intersecting a base-pair position.'''
        problem = im.io.read_npz(im.itu.FAMILY13_STAGE4)
        ibd = problem.info.ibd        
        assert_equal(len(ibd), 60, 'Wrong total number of segments')

        bp = (33554372, 33554373)
        segments = ibd.segments_intersecting_bp_segment(bp)
        assert_equal(len(segments), 28, 'Wrong total number of segments')
        assert_equal(len(set([y[0] for x in segments for y in x.samples])), 15, 'Wrong # samples appearing in IBD segments')
        # This call is supposed to do the same thing
        assert_equal(len(im.imputation.ibd_lookup.samples_ibd_at_bp(ibd, bp)), 15, 'Wrong # samples appearing in IBD segments')

    def test_read_snps_to_impute(self):
        '''Test loading SNPs to be imputed from an NPZ file generated with cgi2plink
        (by way of io_genotype.write()).'''
        a = ImputationSet.load(im.itu.IMPUTE_RARE)
        g = a.genotype
        assert_size_equals(g, 146, 98)

#    def test_impute(self):
#        '''Test imputing SNPs on chromosome 22.'''
#        p = im.io.read_npz(im.itu.IMPUTE_PHASED_SAMPLE)
#        ibd = im.smart_segment_set.SmartSegmentSet.load(1415, im.itu.IMPUTE_IBD_SAMPLE)
#        assert_equal(len(ibd), 1415, 'Wrong total number of entries')
#        assert_equal(ibd.size, 5492, 'Wrong total number of segments')
#        t = im.imputation.ImputationSet.from_file(p.pedigree, im.itu.IMPUTE_RARE_SAMPLE) #@UndefinedVariable
#        snps = np.where(t.snp['chrom'] == 22)[0] # SNP list out of all SNPs in t to impute
#        im.imputation.iibd.impute(p.haplotype, ibd, t, snp=snps, debug=0)
#        t_imputed = t.imputed_data[snps, :, :]
#        assert_almost_equal((1.0 * len(np.where(t_imputed)[0])) / t_imputed.size, 0.3125, decimal=3, err_msg='Wrong imputation call rate')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
