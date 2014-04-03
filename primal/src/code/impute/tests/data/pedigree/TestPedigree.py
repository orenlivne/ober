'''
============================================================
Test class Pedigree basic operations. 

Created on May 30, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''

import networkx as nx
import unittest
from impute.data.Pedigree import Pedigree
from numpy.testing.utils import assert_equal
from impute import impute_test_util as itu

class TestPedigree(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_create_empty(self):
        '''Create an empty pedigree instance.'''
        p = Pedigree(nx.DiGraph())
        assert_equal(p._graph.number_of_nodes(), 0, 'Incorrect number of nodes')

    def test_load_hut(self):
        '''Load hutterites pedigree from file.'''
        # Load data from text file to compare with the load result
        p = itu.read_pedigree_from_test_file(itu.HUTT_PED)
        assert_equal(p.graph.number_of_edges(), 7200, 'Incorrect number of edges')
        assert_equal(p.n, 3671, 'Incorrect number of nodes')

    def test_load_hut_and_recode(self):
        '''Load hutterites pedigree from file; recode to consecutive IDs based on a genotyped
        ID list.'''
        p = itu.read_pedigree_from_test_file(itu.HUTT_PED, genotyped_id_file=itu.GENOTYPE_SAMPLE + '.tfam')
        assert_equal(p.graph.number_of_edges(), 7200, 'Incorrect number of edges')
        assert_equal(p.n, 3671, 'Incorrect number of nodes')
        assert_equal(p.graph.nodes(), range(0,p.n), 'Incorrect node IDs')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
