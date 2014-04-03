'''
============================================================
General-purpose testing utilities.

Created on May 31, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, time
import numpy as np
from unittest.suite import TestSuite
from unittest.loader import TestLoader
from nose.tools import nottest
from numpy.ma.testutils import assert_equal, assert_almost_equal
from timeit import Timer

#---------------------------------------------
# Constants
#---------------------------------------------
TEST_DIR = os.environ['TEST_DATA_DIR']

#---------------------------------------------
# Methods
#---------------------------------------------

@nottest
def abs_path(file_name):
    '''Return the fully-qualified path of the relative test file name file_name.'''    
    return TEST_DIR + '/' + file_name

@nottest
def load_tests_from_classes(test_cases):
    suite = TestSuite()
    loader = TestLoader()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite

@nottest
def timer(title, f, verbose=True, *args, **kwargs):
    '''Time a function call with a return value.'''
    start = time.time()
    result = f(*args, **kwargs)
    if verbose:
        print '%s time: %.3f sec' % (title, time.time() - start)
    return result        

@nottest
def timer_no_retval(title, f, verbose=True):
    '''Time a function call without a return value.
    WARNING: seems fishy. Is f really called?'''
    start = time.time()
    f()
    if verbose:
        print '%s time: %.3e sec' % (title, time.time() - start)

@nottest
def time_of_call(stmt, times=1):
    '''Time running a function(args) <times> times. Returns the average time per call.'''
    return Timer(stmt).timeit(times)/(1.0*times)

#---------------------------------------------
# Assertions
#---------------------------------------------

@nottest
def assert_equal_as_sets(actual, desired, err_msg=''):
    '''Assert that two lists are equal up to a permutation.'''
    assert_equal(set(actual), set(desired), err_msg=err_msg)

@nottest
def assert_positions_almost_equal(actual, expected, decimal=10, err_msg='', allow_flip=False):
    '''Assert that node drawing positions are correct up to symmetry across the y-axis
    (if allow_flip=True).'''
    k = expected.iterkeys().next() 
    if allow_flip and actual[k][0] * expected[k][0] < 0:
        actual = dict((k,(-x,y)) for (k,(x,y)) in actual.iteritems())
    # Maintain the same key order in both dictionaries when comparing
    keys = list(expected.keys()) 
# Debugging printouts
#    a = np.concatenate((np.array([np.array(keys)]).transpose(), 
#                          np.array([actual[k] for k in keys]),
#                          np.array([expected[k] for k in keys])),  
#                          axis=1)
#    for i in a:
#        print '%5d %.10f %.10f %.10f %.10f' % (i[0], i[1], i[2], i[3], i[4])
    assert_almost_equal(np.array([actual[k] for k in keys]), np.array([expected[k] for k in keys]), 
                        decimal=decimal, err_msg=err_msg)
