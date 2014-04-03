'''
============================================================
Python class syntax unit tests. 

Created on May 30, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
import networkx as nx
import unittest
from tests.SimpleClasses import Angle, Parrot, Clazz1, Clazz2, Seq, Slicer
from numpy.ma.testutils import assert_equal

class TestBasicClass(unittest.TestCase):
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_simple_assertion(self):
        s = 'hello, world'
        assert len(s) == 12, 'Incorrect string size'

    def test_graph_creation(self):
        n = 10
        g = nx.path_graph(n, nx.DiGraph())
        assert len(g.nodes()) == n, 'Incorrect number of nodes'

    def test_property(self):
        angle = Angle(10)
        assert angle.rad == 10, 'Property not initialized properly'
        angle.rad = 20
        assert angle.rad == 20, 'Property not initialized properly'
        del angle.rad
        try:
            assert angle.rad == 20, 'Property not deleted properly'
        except:
            pass
        
        c = Clazz2()
        assert c.x == 10, 'Property not initialized properly'
        
        p = Parrot()
        assert p.voltage == 100000, 'Property not initialized properly'
        
        c = Clazz1()
        #print 'voltage', c.voltage, 'x', c.x
        assert c.voltage == 100000, 'Property not set properly'
        assert c.x == 200, 'Property not set properly'
        c.x = 10
        #print 'voltage', c.voltage, 'x', c.x
        assert c.voltage == 100000, 'Property not set properly'
        assert c.x == 10, 'Property not set properly'
        
    def test_iterator(self):
        '''Test an iterator class.'''
        s = Seq()
        n = 0
        for i in s:
            n += 1
            if n == 10:
                break
        assert i == 100, 'Wrong iteration result'

    def test_slicer(self):
        '''Test slice arguments for indexing.'''
        a = np.array([[[1,2,3], [4,5,6]],[[7,8,9], [4,5,6]]])
        s = Slicer(a, 0)
        assert_equal(s[:,:], a[0,:,:], 'Wrong slicing result')
        assert_equal(s[:,2], a[0,:,2], 'Wrong slicing result')
        assert_equal(s[1,2], a[0,1,2], 'Wrong slicing result')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
