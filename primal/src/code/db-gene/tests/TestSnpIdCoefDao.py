'''
============================================================
Test retrieving Hutterities kinship information. 

Created on December 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np
from numpy.ma.testutils import assert_equal, assert_almost_equal
from sqlalchemy import create_engine
from db_gene.snp.snp_db_dao import DEFAULT_SNP_DB_DAOS
from db_gene import DEFAULT_URL

class TestSnpIdCoefDao(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Use a localhost UCSC copy.'''
        self.engine = create_engine(DEFAULT_URL)
        self.idcoef_dao = DEFAULT_SNP_DB_DAOS().idcoef_dao
        # Base.metadata.create_all(self.engine) 
        
    def tearDown(self):
        '''Drop the database.'''
        pass
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_num_records(self):
        '''Test getting the total # of SNP records.'''
        assert_equal(self.idcoef_dao.num_records(), 1415 * 1416 / 2, 'Wrong row count from SNP table')

    def test_delta(self):
        '''Test retrieving condensed identity coefficient for a pair of samples.'''
        lam, delta = self.idcoef_dao.id_coefs(168042, 163352)
        assert_almost_equal(lam, 0.341689,
                            decimal=5, err_msg='Wrong SNP lambda coefficient loading from db')
        assert_almost_equal(delta, 
                            np.array([6.89670924e-05,  0.0016811301,    0.00197529634,   0.0418526094,
                                      0.00258931545,   0.0338666878,    0.00123981348,   0.0908498468,
                                      0.825876333]),
                            decimal=5, err_msg='Wrong SNP Delta loading from db')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
