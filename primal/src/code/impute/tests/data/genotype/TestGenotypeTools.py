'''
============================================================
Test basic genotype data structures and low-level algorithms. 

Created on July 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, itertools as it
from impute import impute_test_util as itu
from impute.tools import pedigree_tools as pt
from impute.tools import genotype_tools as gt
from util import is_member, has_at_least_n_members
from numpy.ma.testutils import assert_equal, assert_almost_equal
from test_util import time_of_call

class TestGenotypeTools(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    @classmethod
    def setUpClass(cls):
        '''Load hutterites data only once per the entire test suite.'''
        cls.problem = itu.Templates.problem_hut()
        cls.pedigree = cls.problem.pedigree
        cls.p = cls.pedigree.graph
        cls.genotype = cls.problem.genotype
        cls.g = cls.genotype.sample_id

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_is_member(self):
        '''Check simple set membership logic.'''
        self.assertTrue(is_member(self.g, [self.pedigree.node_of[126251]]),
                        'Should have been in genotyped list but is not')
        self.assertFalse(is_member(self.g, [-1]), 'Should have not been in genotyped list but is')

        self.assertTrue(has_at_least_n_members(set([1, 2, 3]), [1, 2], 1))
        self.assertTrue(has_at_least_n_members(set([1, 2, 3]), [1, 2], 2))
        self.assertTrue(has_at_least_n_members(set([1, 2, 3]), [1, 2], 3))

    def test_in_extended_family(self):
        '''
        Check the logic of the nulcear family membership algorithm on a sample of
        genotyped hutterites and their pedigree. Does not treat polygamous families correctly.
        '''
        # print len(self.g), self.g.__class__, self.g
        assert_almost_equal((1.0 * sum(self.p.in_degree().itervalues())) / self.p.number_of_nodes(),
                            1.961, 3, 'Wrong average pedigree degree')
        self.assertTrue(TestGenotypeTools.__in_extended_family(self.p, self.g, self.pedigree.node_of[106592]),
                        'Should have been in nuclear family')
        assert_equal(self.genotype.num_samples,
                         1415, 'Unexpected total # of genotyped persons')
        assert_equal(sum(TestGenotypeTools.__in_extended_family(self.p, self.g, x) for x in self.g),
                         736, 'Unexpected # of nuclear family members')
        assert_equal([sum(TestGenotypeTools.__in_extended_family(self.p, self.g, x, n) for x in self.g) 
                      for n in np.arange(13, 0, -1)],
                     [736, 736, 736, 752, 780, 837, 870, 934, 977, 1027, 1067, 1122, 1149],
                         'Unexpected # of nuclear family members')
    
    def test_complete_haplotype_full(self):
        '''Test completing a haplotype with one known entry using a fully-called genotype.'''
        h = np.array([[1, 0],
                      [1, 0],
                      [2, 0],
                      [2, 0],
                      [0, 1],
                      [0, 1],
                      [0, 2],
                      [0, 2],
                      [0, 1]])
        g = np.array([[1, 2],
                      [2, 1],
                      [1, 2],
                      [2, 1],
                      [1, 2],
                      [2, 1],
                      [1, 2],
                      [2, 1],
                      [0, 0]])
        h_expected = np.array([[1, 2],
                      [1, 2],
                      [2, 1],
                      [2, 1],
                      [2, 1],
                      [2, 1],
                      [1, 2],
                      [1, 2],
                      [0, 1]])
        h_temp = h.copy()
        gt.complete_haplotype(h_temp, g, np.zeros((g.shape[0],), dtype=int))
        assert_equal(h_temp, h_expected, 'Unexpected haplotype completion')
        
    def test_complete_genotype_partial(self):
        '''Test completing a haplotype with one known entry vs. a partially-called genotype.'''
        self._test_complete_genotype_partial([1, 0], [[1, 0], [1, 0], [1, 2],
                                                      [1, 0], [1, 1], [1, 2],
                                                      [1, 2], [1, 2], [1, 0]])
        self._test_complete_genotype_partial([0, 1], [[0, 1], [0, 1], [2, 1],
                                                      [0, 1], [1, 1], [2, 1],
                                                      [2, 1], [2, 1], [0, 1]])
        self._test_complete_genotype_partial([2, 0], [[2, 0], [2, 1], [2, 0],
                                                      [2, 1], [2, 0], [2, 1],
                                                      [2, 0], [2, 1], [2, 2]])
        self._test_complete_genotype_partial([0, 2], [[0, 2], [1, 2], [0, 2],
                                                      [1, 2], [0, 2], [1, 2],
                                                      [0, 2], [1, 2], [2, 2]])
        
    def test_genotyped_trios(self):
        '''Find all genotyped trios.'''
        trios = self.problem.trios(genotyped=False)
        assert_equal(trios.shape, [3600, 3], 'Unexpected # of ungenotyped trios')
        trios = self.problem.trios(genotyped=True)
        assert_equal(trios.shape, [869, 3], 'Unexpected # of genotyped trios')
    
        f = self.problem.find_families_by_mother(963, genotyped=True)[0]
        assert_equal(f.num_children, 1, 'Unexpected # genotyped family children')
        f = self.problem.find_families_by_mother(963, genotyped=False)[0]
        assert_equal(f.num_children, 3, 'Unexpected # ungenotyped family children')

    def test_in_family_via_trios(self):
        '''A more exact computation of nuclear family membership via genotyped trios. This treats
        polygamous families correctly, unlike in_family(). This version calculates membership
        so that there are at least min_children children (no relation to size of pedigree family).'''
        assert_equal([len(self.problem.families_union(min_children=n)) for n in np.arange(14, 0, -1)],
                     [0, 15, 15, 15, 63, 106, 233, 320, 530, 692, 860, 978, 1075, 1151],
                     'Unexpected # of nuclear family members')
        
    def test_genotyped_trios_caching(self):
        '''Find all genotyped trios and make sure caching them is fast.'''
        def f():
            assert_equal(self.problem.trios().shape, [869, 3], 'Unexpected # of genotyped trios')

        t1 = time_of_call(f, 1)
        # Since the genotyped trios is cached, the second call should be much faster 
        t2 = time_of_call(f, 1)
        self.assertTrue(t2 < 10 * t1, 'Expecting much faster access upon property caching')

    def test_family_dataset(self):
        '''Test the size and number of nuclear families in the single nuclear family data set.'''
        # print len(self.g), self.g.__class__, self.g
        problem = itu.Templates.problem_family(itu.FAMILY7)
        min_children = 3
        itu.assert_size_equals(problem.genotype, 3218, 9)
        assert_equal(len(problem.families_union(min_children=min_children)), 9)
        family = problem.families(min_children)[0]
        assert_equal(family.father, 0, 'Wrong mother ID') 
        assert_equal(family.mother, 1, 'Wrong mother ID')
        assert_equal(family.children, set([2, 3, 4, 5, 6, 7, 8]), 'Wrong children set')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def _test_complete_genotype_partial(self, h, h_expected):
        '''Test completing a haplotype with one known entry vs. a partially-called genotype. Comprehensive
        checks for a single h and all possible genotypes g.'''
        g = np.array(list(it.product(xrange(3), xrange(3))))
        h_temp = np.tile(h, (9, 1)).copy()
        gt.complete_haplotype_partial(h_temp, g)
        assert_equal(h_temp, h_expected, 'Unexpected haplotype completion')
    
    # Method previously in genotyped_tools, replaced by genotyped_trios and family calls
    @staticmethod
    def __in_extended_family(p, g, x, min_children=None):
        '''Check whether a member x of a genotyped-person list g in is in a nuclear family or not according
        to the pedigree graph p.
        This is defined as:
        CASE 1: x's parents are in G and all their children are in G;
        OR 
        CASE 2: x has children, they are all in G and all their parents are in G.
        
        If min_children is specified (not None), the following alternative definition is used instead:
        CASE 1: x's parents are in G and at least min_children of their children are in G;
        OR 
        CASE 2: x has children, of which at least min_children are in G; all of x's children's parents are in G.
        
        Note: this code does not treat polygamous families correctly and assumes that such families are
        'extended', i.e., have more than 2 parents.
        '''
        
        # Optimization: don't bother checking further if x is not in G
        if not p.has_node(x):
            raise Exception('Node ' + str(x) + ' not in pedigree')  
        if not is_member(g, [x]):
            return False
        
        # Case 1
        parents = p.predecessors(x)
        if (is_member(g, parents) and 
            has_at_least_n_members(g, pt.all_successors(p, parents), min_children)):
            return True
        
        # Case 2
        children = p.successors(x)
        if (children and has_at_least_n_members(g, children, min_children) 
            and is_member(g, pt.all_predecessors(p, children))):
            return True
        
        return False
