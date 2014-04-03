'''
============================================================
Impute project - main test suite.

Created on Jun 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import test_util as tu
from tests.TestFamPlot import TestFamPlot

def suite():
    return tu.load_tests_from_classes([TestFamPlot])
