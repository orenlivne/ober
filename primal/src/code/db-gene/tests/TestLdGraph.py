'''
============================================================
Test LD graph algorithms. 

Created on November 21, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, networkx as nx, numpy as np
from numpy.ma.testutils import assert_equal
from db_gene.snp import ld_graph
import util

class TestLdGraph(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_greedy_coloring(self):
        '''Test the greedy coloring algorithm on simple chain and near-chain graphs.'''
        n = 10
        g = nx.path_graph(n)
        
        # Increasing positions vs. node number
        self.__test_greedy_coloring(g, np.arange(n) * 2, [0, 1, 0, 1, 0, 1, 0, 1, 0, 1])
        
        # Decreasing positions vs. node number
        self.__test_greedy_coloring(g, -np.arange(n) * 2, [1, 0, 1, 0, 1, 0, 1, 0, 1, 0])
        
    def test_greedy_coloring_nontrivial_graph(self):
        '''Test the greedy coloring algorithm on simple chain and near-chain graphs.'''
        n = 10
        g = nx.path_graph(n)
        # Make graph non-chain
        g.add_edge(0, 2)
        g.add_edge(0, 3)
        self.__test_greedy_coloring(g, np.arange(n) * 2, [0, 1, 2, 1, 0, 1, 0, 1, 0, 1])
    
    def test_katie_graph(self):
        '''Test an LD graph from Katie that seems to be probelmatic.'''
        g = nx.Graph()
        g.add_nodes_from(['70', '60', '75', '67', '77', '63', '76', '89', '37'])
        g.add_edges_from([('70', '63'), ('70', '60'), ('70', '75'), ('70', '67'), ('70', '77'), ('70', '76'), ('70', '89'), ('70', '37'), ('60', '89'), ('60', '75'), ('60', '67'), ('60', '76'), ('60', '63'), ('75', '67'), ('75', '77'), ('75', '76'), ('75', '89'), ('75', '63'), ('75', '37'), ('67', '77'), ('67', '76'), ('67', '89'), ('67', '63'), ('67', '37'), ('77', '89'), ('77', '76'), ('77', '63'), ('77', '37'), ('63', '76'), ('63', '89'), ('63', '37'), ('76', '89'), ('76', '37')])
        position = dict([('37', 37441365), ('60', 37441980), ('63', 37482151), ('67', 37487632), ('70', 37487723), ('75', 37487866), ('76', 37487873), ('77', 37487975), ('89', 37488499)])

        expected_colors = [0, 0, 1, 2, 3, 4, 5, 6, 7]
        c = ld_graph.greedy_coloring(g, position)
        assert_equal(list(bad_edges(g, c)), [], 'Edges found within a frame')
        assert_equal(map(c.get, sorted(c)), expected_colors, 'Wrong greedy coloring')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __test_greedy_coloring(self, g, position, expected_colors):
        '''Test the greedy coloring algorithm.'''
        c = ld_graph.greedy_coloring(g, dict(zip(g.nodes(), position)))
        assert_equal(list(bad_edges(g, c)), [], 'Edges found within a frame')
        assert_equal(map(c.get, sorted(c)), expected_colors, 'Wrong greedy coloring')

def group_by_value(d):
    '''Convert a dictionary to a dictionary where each key is a d-value and the value is the list of keys
    in d that have this value.'''
    return util.to_set_dict((v, k) for k, v in d.iteritems())

def bad_edges(g, c):
    '''Check that there are no ''bad'' edges within a frame.'''
    frames = group_by_value(c)
    for frame in frames.itervalues():
        for x in frame:
            for y in frame:
                if g.has_edge(x, y): yield x, y
