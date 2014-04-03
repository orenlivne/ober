'''
============================================================
Test the filter chain framework. 

Created on July 5, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''

import unittest
from numpy.ma.testutils import assert_equal
from chain import Filter, FilterChain
from util import Struct

class TestChain(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
 
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_simple_chain(self):
        '''Check that a chain with printout filters works.'''
        chain = FilterChain([IncrementFilter(1), IncrementFilter(2), IncrementFilter(3)])
        request = Struct(x=0)
        self.assertFalse(chain.handle(request))
        assert_equal(request.x, 6, 'Bad chain processing')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------

####################################################################################
class IncrementFilter(Filter):
    '''A filter that adds a value to its a request with the integer field 'x'.'''
    def __init__(self, increment, next_filter=None):
        '''Initialize a filter.'''
        Filter.__init__(self, next_filter)
        self.increment = increment
     
    def handle(self, request):
        '''Print some info.'''
        request.x = request.x + self.increment 
