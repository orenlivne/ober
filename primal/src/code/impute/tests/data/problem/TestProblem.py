'''
============================================================
Test class Problem's basic methods. 

Created on August 13, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
import unittest
import test_util as tu
from impute import impute_test_util as itu
from numpy.testing.utils import assert_equal
from numpy.ma.testutils import assert_not_equal
from impute.data import io

class TestProblem(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_select_samples(self):
        '''Create sub-problem by selecting samples/family..'''
        p = io.read_npz(itu.FAMILY13_STAGE3)
        f = p.families()[0]
        ps = p.sub_problem(f.member_list)
        pf = p.sub_problem_of_parents(f.father, f.mother)

        # pf has family_info entry, ps doesn't
        assert_equal(ps.pedigree.graph.nodes(), pf.pedigree.graph.nodes(), 'Incorrect pedigree graph nodes')
        assert_equal(set(ps.pedigree.graph.edges()), set(pf.pedigree.graph.edges()), 'Incorrect pedigree graph edges')
        assert_equal(ps.pedigree.sample_id.tolist(), pf.pedigree.sample_id.tolist(), 'Incorrect pedigree sample ID array')
        assert_equal(ps.pedigree, pf.pedigree, 'Incorrect pedigree')
        assert_equal(ps.genotype, pf.genotype, 'Incorrect genotype')
        assert_equal(ps.haplotype, pf.haplotype, 'Incorrect haplotype')
        self.assertTrue(ps.info.family_info != pf.info.family_info, 'Incorrect info family')
        assert_not_equal(ps.info, pf.info, 'Incorrect info')
        assert_not_equal(ps, pf, 'pf should have a family_info entry, ps shouldn''t')

    def test_select_family(self):
        '''Another test of creating a sub-problem by selecting a family..'''
        p = io.read_npz(itu.FAMILY13_STAGE3)
        p.sub_problem_of_parents(0, 1)

    def test_neighbors(self):
        '''Test retrieving node's pedigree neighbors.'''
        p = itu.Templates.problem_hut()
        node = 1298
        
        expected = np.genfromtxt(tu.abs_path('ibd_distant/nbhrs_' + repr(node) + '.txt'), dtype=int).tolist()
        tu.assert_equal_as_sets(p.pedigree.neighbors(node, 4, False), expected, 'Incorrect genotyped neighbor list')
        
        expected = np.genfromtxt(tu.abs_path('ibd_distant/nbhrs_' + repr(node) + '_genotyped.txt'), dtype=int).tolist()
        tu.assert_equal_as_sets(p.pedigree.neighbors(node, 4, True), expected, 'Incorrect genotyped neighbor list')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    