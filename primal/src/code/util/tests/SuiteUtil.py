'''
============================================================
Misc code project - main test suite.

Created on Jun 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import test_util as tu
from tests.TestChain import TestChain
from tests.TestBasicOperations import TestBasicOperations
from tests.TestBasicClass import TestBasicClass
from tests.TestBinarySearchTree import TestBinarySearchTree
from tests.TestSavez import TestSavez
from tests.TestPil import TestPil
from tests.TestHashable import TestHashable
from tests.TestItemUtil import TestItemUtil
from tests.TestStatUtil import TestStatUtil
from tests.db.TestSqlite3 import TestSqlite3
from tests.db.TestSqlAlchemy import TestSqlAlchemy
from tests.math.TestIndexUtil import TestIndexUtil
from tests.TestUtilNeighbor import TestUtilNeighbor
from tests.hmm.TestHmmt import TestHmmt
from tests.hmm.TestHmm import TestHmm
from tests.TestIntervalTree import TestIntervalTree
from tests.TestSingleton import TestSingleton
from tests import TestRobust

def suite():
    return tu.load_tests_from_classes([TestChain, TestBasicOperations, TestBasicClass,
                                       TestBinarySearchTree, TestSavez,
                                       TestPil, TestHashable, TestItemUtil,
                                       TestStatUtil,
                                       TestSqlite3, TestSqlAlchemy,
                                       TestHmm, TestHmmt,
                                       TestIndexUtil, TestUtilNeighbor, TestIntervalTree,
                                       TestSingleton,
                                       TestRobust])
