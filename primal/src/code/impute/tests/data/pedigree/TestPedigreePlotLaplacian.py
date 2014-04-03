'''
============================================================
Test class PedigreeTools - pedigree algorithms. 

Created on May 30, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, test_util, numpy as np, networkx as nx
from impute.tools import pedigree_plot_laplacian
from nose.tools import nottest
from impute.data import io_pedigree
from impute import impute_test_util as itu
from numpy.ma.testutils import assert_equal
from impute.tools.laplacian import laplacian, sym_adjacency

class TestPedigreePlotLaplacian(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_layout_position(self):
        '''Test family computation.'''
        positions = pedigree_plot_laplacian._layout_positions(nx.path_graph(10, create_using=nx.DiGraph()))
        expected = {0: (-0.44170765403093781 , 9), 
                    1: (-0.39847023129619913 , 8), 
                    2: (-0.31622776601683716 , 7), 
                    3: (-0.20303072371134445 , 6), 
                    4: (-0.069959619570752599, 5), 
                    5: ( 0.069959619570755624, 4), 
                    6: ( 0.20303072371134645 , 3), 
                    7: ( 0.3162277660168385  , 2), 
                    8: ( 0.39847023129620007 , 1), 
                    9: ( 0.44170765403093859 , 0)}
        test_util.assert_positions_almost_equal(positions, expected, decimal=10, 
                                                err_msg='Wrong positions', allow_flip=True)

    def test_marriage_graph(self):
        '''Test generating the extended graph.'''
        p = io_pedigree.read(itu.HUTT_PED)
        g = p.graph
        g_extended = pedigree_plot_laplacian._marriage_graph(g)
        
        assert_equal(g.number_of_nodes(), 3671, 'Wrong number of nodes')
        assert_equal(g.number_of_edges(), 7200, 'Wrong number of edges')
        
        assert_equal(g_extended.number_of_nodes(), 4661, 'Wrong number of nodes')
        assert_equal(g_extended.number_of_edges(), 5580, 'Wrong number of edges')
        
    def test_marriage_graph_layout_positions(self):
        '''Test generating the extended graph layout position.'''
        p = io_pedigree.read(itu.SMALL_FILE)
        g = p.graph
        g_extended = pedigree_plot_laplacian._marriage_graph(g)
        
        assert_equal(g.number_of_nodes(), 8, 'Wrong number of nodes')
        assert_equal(g.number_of_edges(), 10, 'Wrong number of edges')
        
        assert_equal(g_extended.number_of_nodes(), 12, 'Wrong number of nodes')
        assert_equal(g_extended.number_of_edges(), 13, 'Wrong number of edges')

        #positions = 
        pedigree_plot_laplacian._layout_positions(g, g_extended)
        #expected = 
        {1: (-0.23562670914672229, 2), 2: (-0.063382627268225591, 3), 3: (-0.23562670914672237, 3), 4: (0.1502736499569044, 1), 5: (0.15027364995690434, 2), 6: (0.43352532974526942, 0), 7: (0.43352532974526858, 0), 8: (-0.48586913569278134, 2), -1: (-0.054615070741569148, 2.5), -4: (0.31358843226722749, 0.5), -3: (-0.054615070741569086, 2.5), -2: (-0.35145106893398231, 1.5)}
        # Positions seem to  slightly vary from one platform to another 
        #test_util.assert_positions_almost_equal(positions, expected, decimal=9, 
        #                                        err_msg='Wrong positions', allow_flip = True)

    @nottest
    def too_slow_test_layout_hut(self):
        '''Test generating layout coordinates for the Hutterites pedigree. This test is more about
        speed, to see that we do it in reasonable time.'''
        p = io_pedigree.read(itu.HUTT_PED)
        positions = pedigree_plot_laplacian._layout_positions(p.graph) #@UnusedVariable
        # Save to file - slow
        '''
        print positions
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        nx.draw_networkx(p.graph, pos=positions)
        plt.savefig('hutterites.png')
        '''
        
    def test_matrices(self):
        '''Test adjacency and Laplacian matrix computation.'''
        # Load data from text file to compare with the load result
        p = itu.read_pedigree_from_test_file(itu.HUTT_PED)
        A = sym_adjacency(p.graph)
        L = laplacian(p.graph)
        assert not np.nonzero(A-A.T)[0], 'Adjacency matrix returned from sym_adjacency() should be symmetric'
        assert not np.nonzero(L-L.T)[0], 'Laplacian matrix should be symmetric'
        assert (abs(L.sum(axis=1)) < 1e-8).all(), 'Laplacian matrix should be zero-row-sum'
        #lam, v = eigsh(L, 'SM', 1)
