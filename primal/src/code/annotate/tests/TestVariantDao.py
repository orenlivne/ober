'''
============================================================
Test retrieving Gene information from a local mirror of the
UCSC gene browser. 

Created on November 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, StringIO
from numpy.ma.testutils import assert_equal
from annotate.hutt_dao import VariantDao, REGION_CATEGORY_OF,\
    VariantSummaryReport
from collections import Counter
# from sqlalchemy import create_engine

class TestHuttDao(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    TEST_URL = 'mysql://hutt:hutt@127.0.0.1/hutt'  # Should be built from $OBER/testdata/annotate/data.mysql.qc.10000 using the create-annotation-db script
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Use a localhost UCSC copy.'''
        self.variant_dao = VariantDao(TestHuttDao.TEST_URL)
        
    def tearDown(self):
        '''Drop the database.'''
        pass
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_num_records(self):
        '''Test getting the total # of SNP records.'''
        assert_equal(self.variant_dao.num_records(), 10000, 'Wrong row count from SNP table')

    def test_variant_count(self):
        '''Test variant summary counts queries.'''
        maf = np.concatenate((np.arange(0, 0.11, 0.01), [0.5]))
        for i in xrange(len(maf) - 1):
            result = Counter(dict((REGION_CATEGORY_OF[k], v) for k, v in self.variant_dao.variant_count_report_maf_bin('region', maf[i], maf[i + 1])))
            self.assertTrue(max(result.itervalues()) > 0, 'No non-trivial variant category count found in result dictionary')
                
    def test_variant_count_report_region(self):
        '''Test variant summary counts queries, grouped by region.'''
        self._test_variant_count_report('region')
                        
    def test_variant_count_report_coding(self):
        '''Test variant summary counts queries, grouped by coding.'''
        self._test_variant_count_report('coding')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def _test_variant_count_report(self, group):
        '''Test variant summary counts queries, grouped by region.'''
        report = self.variant_dao.variant_count_report(group)
        output = StringIO.StringIO()
        report.save(output)
        str1 = output.getvalue()
        output.close()

        output = StringIO.StringIO()
        report.save(output)
        report2 = VariantSummaryReport.load(iter(output.getvalue().split('\n')[:-1]))
        output.close()
        output = StringIO.StringIO()
        report2.save(output)
        str2 = output.getvalue()
        assert_equal(str1, str2, 'Saving-loading-saving doesn\'t give the same serialization as saving')
