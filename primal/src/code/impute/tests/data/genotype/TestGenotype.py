'''
============================================================
Test class Pedigree basic operations. 

Created on May 30, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''

import numpy as np
import unittest
from impute import impute_test_util as itu
from numpy.testing.utils import assert_equal
from impute.data.factory import GenotypeFactory
from impute.data import io_genotype

class TestGenotype(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_create_from_mock_data(self):
        '''Create a simple genotype set from the hutterites pedigree and some mock genotype data.'''
        # Load data from text file to compare with the load result
        snp = np.array([(0, 'rs1', 0., 12), 
                        (0, 'rs2', 0., 34), 
                        (0, 'rs3', 0., 56),
                        (0, 'rs4', 0., 78)],
                      dtype={'names': ('chrom', 'snp', 'dist_cm', 'base_pair'), 
                             'formats': ('i2', 'S12', 'i8', 'i8')})               
        sample_id = [126251, 111161]
        data = np.array([[[1, 2]], [[2, 2]], [[1, 2]], [[1, 1]]])
        g = GenotypeFactory.new_instance('genotype', data, snp, sample_id)
        itu.assert_size_equals(g, 4, 1)
        assert_equal(4, g.num_snps, 'Incorrect number of SNPS')
        assert_equal(g.segment_intersect([0, 40]), [0, 2], 'Wrong interval intersection')
        assert_equal([0, 2], g.segment_intersect([10,40]), 'Wrong interval intersection')
        assert_equal([0, 3], g.segment_intersect([10,60]), 'Wrong interval intersection')
        assert_equal([1, 3], g.segment_intersect([20,60]), 'Wrong interval intersection')
        assert_equal([0, 4], g.segment_intersect([0,100]), 'Wrong interval intersection')
        assert_equal([1, 4], g.segment_intersect([20,100]), 'Wrong interval intersection')

    def test_create_from_file(self):
        g = io_genotype.read('plink', 'genotype', prefix=itu.GENOTYPE_SAMPLE)
        self.__assert_genotype_set_stats_correct(g)
        self.__assert_snp_equals((22, 'rs1654', 0, 17596388), g.snp[7])
                
    def test_create_from_file_no_sample_ids(self):
        g = io_genotype.read('plink', 'genotype', tped=itu.GENOTYPE_SAMPLE+'.tped', load_ids=False)
        self.__assert_genotype_set_stats_correct(g)


    def test_nearest_snp(self):
        '''Test finding the nearest SNP of a given base pair location.'''
        g = io_genotype.read('plink', 'genotype', prefix=itu.GENOTYPE_SAMPLE)
        assert_equal(list(g.nearest_snp_multiple(g.base_pair)), g.snp_range, 'Wrong nearest SNP location')
        i = 5
        assert_equal(g.nearest_snp(g.base_pair[i]+0.01), i, 'Wrong nearest SNP location')
        assert_equal(g.nearest_snp(g.base_pair[i]-0.01), i, 'Wrong nearest SNP location')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------

    def __assert_genotype_set_stats_correct(self, g):
        (n, snps) = (1415, 8)
        itu.assert_size_equals(g, snps, n)
        assert_equal(len(g.sample_id), n, 'Incorrect sample ID set size')
        assert_equal(g.num_snps, snps, 'Incorrect number of SNPS')
        assert_equal(g.data.shape, [snps, n, 2], 'Incorrect genotype data array shape')

    def __assert_snp_equals(self, actual, snp):
        assert_equal(actual[0], snp['chrom'], 'Wrong SNP metadata')
        assert_equal(actual[1], snp['name'], 'Wrong SNP metadata')
        assert_equal(actual[2], snp['dist_cm'], 'Wrong SNP metadata')
        assert_equal(actual[3], snp['base_pair'], 'Wrong SNP metadata')