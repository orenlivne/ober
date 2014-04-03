'''
============================================================
Test numpy binary serialization. 

Created on Aug 4, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''

import unittest
import networkx as nx
from numpy.ma.testutils import assert_equal
import numpy as np
import tempfile

class TestSavez(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
 
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_objects(self):
        '''Check that saving and loading using numpy binary format works. As long as the
        object is wrapped in a numpy array, this works.'''
        self.__test_save_and_load(np.array([1,2,3]))
        self.__test_save_and_load(np.array(['str']))
        self.__test_save_and_load(np.array(ComplexObject('x',y=np.array([1,2,3]))))
        
    def test_graph(self):
        self.__test_save_and_load(np.array([nx.path_graph(10, nx.DiGraph())]))
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __test_save_and_load(self, x):
        out_file = tempfile.TemporaryFile()
        np.savez(out_file, x=x)
        out_file.seek(0) # Only needed here to simulate closing & reopening file
        x2 = np.load(out_file)
        y = x2['x']
        assert_equal(x, y, 'Saving and loading did not restore the original object')
        #out_file.close() # Would be needed with a normal file descriptor

####################################################################################
class ComplexObject(object):
    '''A complex object to serialize.'''
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
    def __key(self):
        '''Uniquely-identifying key of this family: parent tuple.'''
        return (self.x, self.y.tolist()) # numpy array hash comparison does not seem to work as expected

    def __eq__(self, other):
        '''Equality of families.'''
        return self.__key() == other.__key()
    
    def __ne__(self, y):
        '''Inequality of objects.'''
        return self.__key() != y.__key()

    def __hash__(self):
        '''Object hash code.'''
        return hash(self.__key())

    def __repr__(self):
        return 'ComplexObject[x=' + repr(self.x) + ', y=' + repr(self.y)+']' 
