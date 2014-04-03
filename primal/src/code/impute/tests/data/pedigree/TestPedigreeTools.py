'''
============================================================
Test class PedigreeTools - pedigree algorithms. 

Created on May 30, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, tempfile, unittest, impute.tools.pedigree_tools as pt
from numpy.testing.utils import assert_equal
from networkx.algorithms.shortest_paths.generic import shortest_path
from numpy.core.numeric import Infinity
from impute import impute_test_util as itu
from impute.data import io_pedigree
from tests.data.pedigree.pedigree_old_study import PedigreeOldStudyReader
from impute.data.constants import PATERNAL, MATERNAL
from impute.data.io import read_npz
from collections import OrderedDict

class TestPedigreeTools(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_check_generation_number(self):
        '''
        Load hutterites pedigree from file including the old_generation numbers
        given in that file. Check if each child's parents have the same old_generation #. Such
        an ordering cannot be guaranteed to exist for any DAG, but maybe here.
        '''
        file_name = itu.HUTT_PED
        p = PedigreeOldStudyReader().read(file_name)

        # Check that generation order from study is consistent
        self.__check_generation_order_consistent(p, p.old_generation)
        self.__assert_max_generation_gap_equals(p, p.old_generation, 7) 
        
        # Test our algorithm vs. old study algorithm without alien node alignment. Only
        # a single node should be different (an alien, 999999981)
        depth = pt.node_depth(p.graph, max_generation_gap=0)
        self.__test_generation_order_differs_by(p, depth, 1)
        self.__assert_max_generation_gap_equals(p, depth, 10) 

        # Test depth calculation with a sufficiently enough max generation gap between parent
        # and child, as is implicitly assumed in Gaixin's and Mark's old study code
        depth = pt.node_depth(p.graph, max_generation_gap=8)
        self.__test_generation_order_differs_by(p, depth, 0)
        self.__assert_max_generation_gap_equals(p, depth, 7)

        # With alien alignment and our default, smaller max generation gap, check that our
        # generation order is consistent, but it should be quite different than the old study's
        depth = pt.node_depth(p.graph)
        self.__test_generation_order_differs_by(p, depth, 19)
        # Our order almost attained the requested (default) max generation depth of 3
        self.__assert_max_generation_gap_equals(p, depth, 4)
        
        '''
        # Print a report of our gen order next to the original in the pedigree file format
        data = np.genfromtxt(file_name, np.dtype(int))
        for i in np.arange(0, data.shape[0]):
            node = data[i,0]
            print '%7d %7d %7d %7d %7d %7d' % (node, data[i,1], data[i,2], data[i,3], depth[node], (depth[node] != data[i,3]))
        '''
        
    def test_families(self):
        '''Test family computation.'''
        p = io_pedigree.read(itu.HUTT_PED)
        assert_equal(len(list(pt.families(p.graph))), 990, 'Wrong number of pedigree families')

    def test_sibs(self):
        '''Find all siblings of a node.'''
        p = itu.Templates.problem_hut().pedigree
        assert_equal(set(pt.sibs(p, 963, PATERNAL)), set([1, 2786, 963, 2788, 973, 1811, 2787]), 'Unexpected sibling set')
        p = itu.Templates.problem_hut().pedigree
        assert_equal(set(pt.sibs(p, 963, MATERNAL)), set([1, 2786, 963, 2788, 973, 1811, 2787]), 'Unexpected sibling set')

    def test_lca_small(self):
        '''Test lowest common ancestor computation in a small pedigree.'''
        p = io_pedigree.read(itu.SMALL_FILE)

        # Direct siblings
        self.__compute_and_check_lca(p, 6, 7, [4, 5], 2)

        # Far siblings    
        u = 6
        v = 8
        w = self.__compute_and_check_lca(p, u, v, [3], 3) 
        assert_equal(shortest_path(p.graph, w, u), [3, 5, 6], 'Wrong shorted path from ancestor to node u')
        assert_equal(shortest_path(p.graph, w, v), [3, 8], 'Wrong shorted path from ancestor to node v')
        
        # No common ancestor exists
        self.__compute_and_check_lca(p, 1, 2, [None], Infinity)

    def test_lca_hut(self):
        '''Test lowest common ancestor computation in a large pedigree.'''
        p = io_pedigree.read(itu.HUTT_PED)
        self.__compute_and_check_lca(p, 169512, 170362, [8551], 9)
    
    def test_draw_pedigree(self):
        '''Test famplot pedigree drawing tool integration.'''
        p = read_npz(itu.NBHRS1298_STAGE4)
        self.__test_draw_pedigree(p)
    
    def test_surrogate_parents(self):
        '''Find all surrogate parents of a node.'''
        p = itu.Templates.problem_hut().pedigree
        g = p.graph
        assert_equal(pt.surrogate_parents(g, 1298, 0), {}, 'Unexpected surrogate parent set')
        assert_equal(pt.surrogate_parents(g, 1298, 1), OrderedDict([(1480, 1), (1719, 1), (1217, 1), (1412, 1), (1413, 1)]), 'Unexpected surrogate parent set')
        assert_equal(pt.surrogate_parents(g, 1298, 2), OrderedDict([(1480, 1), (1719, 1), (2028, 2), (1829, 2), (2064, 2), (2082, 2), (1217, 1), (1412, 1), (1413, 1), (2627, 2), (2628, 2), (2629, 2), (2630, 2), (2755, 2), (2756, 2), (1756, 2)]), 'Unexpected surrogate parent set')
        assert_equal(pt.surrogate_parents(g, 1298, 3), OrderedDict([(1480, 1), (1719, 1), (2028, 2), (1829, 2), (2064, 2), (2082, 2), (2180, 3), (2158, 3), (2155, 3), (2070, 3), (2275, 3), (2278, 3), (2202, 3), (2171, 3), (1217, 1), (1412, 1), (1413, 1), (3236, 3), (2627, 2), (2628, 2), (2629, 2), (2630, 2), (617, 3), (1188, 3), (2755, 2), (1110, 3), (2756, 2), (1756, 2), (2982, 3), (1482, 3), (1484, 3), (1486, 3), (1487, 3), (1491, 3), (1492, 3), (2644, 3), (1728, 3), (134, 3), (1735, 3), (2635, 3), (2770, 3), (1749, 3), (1751, 3), (1752, 3), (1721, 3), (1754, 3), (1724, 3), (1726, 3)]), 'Unexpected surrogate parent set')
        assert_equal(pt.surrogate_parents(g, 1298, 3, min_depth=2), OrderedDict([(2028, 2), (1829, 2), (2064, 2), (2082, 2), (2180, 3), (2158, 3), (2155, 3), (2070, 3), (2275, 3), (2278, 3), (2202, 3), (2171, 3), (3236, 3), (2627, 2), (2628, 2), (2629, 2), (2630, 2), (617, 3), (1188, 3), (2755, 2), (1110, 3), (2756, 2), (1756, 2), (2982, 3), (1482, 3), (1484, 3), (1486, 3), (1487, 3), (1491, 3), (1492, 3), (2644, 3), (1728, 3), (134, 3), (1735, 3), (2635, 3), (2770, 3), (1749, 3), (1751, 3), (1752, 3), (1721, 3), (1754, 3), (1724, 3), (1726, 3)]), 'Unexpected surrogate parent set')

#----------------------------------
    # Add another test with a simpler pedigree to make sure we don't catch spouses
    # Test successors=True
#----------------------------------

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __compute_and_check_lca(self, p, u, v, expected_w_set, expected_d):
        (w, d) = pt.lowest_common_ancestor(p.graph, u, v)
        self.assertIn(w, expected_w_set, 'Wrong LCA: actual ' + repr(w) + ' expected ' + repr(expected_w_set))
        assert_equal(d, expected_d, 'Wrong LCA distance')
        return w

    def __test_generation_order_differs_by(self, p, depth, num_differences):
        '''Test that our generation order depth for the pedigree p is different than reference_depth
        (from old study) at exactly num_differences nodes.'''
        reference_depth = p.old_generation 
        self.__check_generation_order_consistent(p, depth)
        if (num_differences == 0):
            assert depth == reference_depth, 'Our depth should be the same as study''s'
        else:
            assert depth != reference_depth, 'Our depth should not be the same as study''s'
        assert_equal(len(TestPedigreeTools.__depth_difference(depth, reference_depth)),
                     num_differences, '%d depth differences should have been detected' 
                     % num_differences)
    
    def __check_generation_order_consistent(self, p, generation):
        '''Check whether the node list generation is a consistent (topological-sort depth)
        generation number on the pedigree p.
                  
        Check whether the ALL parents are in the previous old_generation of the child.
        Throw an exception if parents do not strictly precede the child, even if they
        are not all in the same old_generation.
        
        This is also a demo of a generator function.
        ''' 
        def check_consistent_order():
            for (child, gen) in generation.iteritems():
                expected = gen - 1
                for parent in p._graph.predecessors_iter(child):
                    parent_gen = generation[parent]
                    if (parent_gen >= gen):
                        msg = 'Parent %d generation %d does not precede child %d generation %d' \
                        % (parent, parent_gen, child, gen)
                        raise ValueError(msg)
                    elif (parent_gen != expected):
                        yield (child, False)
                yield (child, True)
            
        for _ in check_consistent_order():
            pass

    def __assert_max_generation_gap_equals(self, p, depth, expected):
        '''Assert that the maximum generation gap between a parent and a child in the
        generation order depth for the pedigree p equals expected'''
        assert_equal(pt.max_generation_gap(p.graph, depth), expected, 'Wrong max generation gap') 

    def __test_draw_pedigree(self, p):
        '''Test drawing a pedigree as an EPS file using the default node placement.'''
        out = tempfile.NamedTemporaryFile(suffix='.eps', delete=False)
        pt.draw_pedigree(p.pedigree, out.name)
        out.close()
        os.remove(out.name)

    @staticmethod
    def __depth_difference(a, b):
        '''Difference between two depth orders.'''
        return [x for x in a if a[x] != b[x]]
