'''
============================================================
Test retrieving Hutterities kinship information. 

Created on December 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, test_util as tu
from numpy.ma.testutils import assert_equal, assert_almost_equal
from db_gene.snp.file_dao import IdCoefDao, KinshipDao

class TestFileIdCoefDao(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Use a localhost UCSC copy.'''
        self.idcoef_dao = IdCoefDao(tu.abs_path('pedigree/example.id'))
        self.kinship_dao = KinshipDao(tu.abs_path('pedigree/example.kinship'))
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_num_records(self):
        '''Test getting the total # of SNP records.'''
        assert_equal(self.idcoef_dao.num_records(), 3 ** 2, 'Wrong row count from SNP table')

    def test_delta(self):
        '''Test retrieving condensed identity coefficient for a pair of samples.'''
        lam, delta = self.idcoef_dao.id_coefs(168042, 163352)
        assert_almost_equal(lam, 0.341689,
                            decimal=5, err_msg='Wrong SNP lambda coefficient loading from file')
        assert_almost_equal(delta,
                            np.array([6.89670924e-05, 0.0016811301, 0.00197529634, 0.0418526094,
                                      0.00258931545, 0.0338666878, 0.00123981348, 0.0908498468,
                                      0.825876333]),
                            decimal=5, err_msg='Wrong SNP Delta loading from file')

    def test_kinship(self):
        '''Test retrieving kinship coefficients for a pair of samples.'''
        assert_equal(self.kinship_dao.num_records(), 3 ** 2, 'Wrong row count from SNP table')

        _, delta = self.idcoef_dao.id_coefs(168042, 163352)
        f = self.kinship_dao.kinship(168042, 163352)
        
        assert_almost_equal(f, _kinship(delta),
                            decimal=5, err_msg='Kinship coefficient does not match value derived from identity coefficients')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
def _kinship(d):
    '''Calculate the kinship coefficient from identity coefficients.'''
    return d[0] + 0.5 * (d[2] + d[4] + d[6]) + 0.25 * d[7]
