'''
============================================================
Test retrieving Gene information from a local mirror of the
UCSC gene browser. 

Created on November 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest
from numpy.ma.testutils import assert_equal
from sqlalchemy import create_engine
from nose.tools import nottest
from db_gene.ucsc.ucsc_dao import SnpDao

class TestUcscSnpDao(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Use a localhost UCSC copy.'''
        url = 'mysql://ucsc:ucsc@localhost/ucsc'
        self.engine = create_engine(url)
        self.dao = SnpDao(url)
        # Base.metadata.create_all(self.engine) 
        
    def tearDown(self):
        '''Drop the database.'''
        pass
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_num_records(self):
        '''Test getting the total # of SNP records.'''
        assert_equal(self.dao.num_records(), 54212080, 'Wrong row count from SNP table')

    # Slow test; can take ~20 seconds
    @nottest
    def inactive_test_snps_on_chromosome(self):
        '''Test retrieving SNPs by chromosome.'''
        chrom = 22
        assert_equal(len(list(self.dao.chrom_snps(chrom))), 730365, 'Wrong number of SNPs on chromosome %d' % (chrom,))

    def test_get_snps(self):
        '''Test retrieving SNPs by chromosome.'''
        snps = self.dao.get_snps(['rs144773400', 'rs143255646', 'rs145599635'])
        assert_equal(len(snps), 3, 'Wrong number returned SNPs entries')
        assert_equal([snp.bp for snp in snps], [10229, 10145, 10234], 'Wrong SNP base-pair positions')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
