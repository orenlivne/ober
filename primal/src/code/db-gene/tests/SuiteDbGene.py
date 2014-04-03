'''
============================================================
Misc code project - main test suite.

Created on Jun 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import test_util as tu
from tests.TestUcscGeneDao import TestUcscGeneDao
from tests.TestSnpLdDao import TestSnpLdDao
from tests.TestUcscSnpDao import TestUcscSnpDao
from tests.TestLdGraph import TestLdGraph
from tests.TestSnpIdCoefDao import TestSnpIdCoefDao
from tests.TestFileIdCoefDao import TestFileIdCoefDao

def suite():
    return tu.load_tests_from_classes([TestUcscGeneDao, TestUcscSnpDao, TestSnpLdDao,
                                       TestSnpIdCoefDao, TestFileIdCoefDao,
                                       TestLdGraph])
