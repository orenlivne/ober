'''
============================================================
Test phasing within nuclear families: comparing sibs
with no genotyped parents.

Created on July 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest
from impute import impute_test_util as itu
from numpy.ma.testutils import assert_equal
from impute.data import io
from impute.tools.param import PhaseParam
from impute.phasing.phase_distant import family_sib_comparison_phaser

class TestPhaseFamilyCompareSibs(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        # The way to load a pedigree in conjunction with a genotype set is to recode
        # its sample IDs to consecutive for easier access by phasers.

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_family_963(self):
        '''Test comparing sibs with non-genotyped parents (stage 4). This was a problematic family.'''
        problem = io.read_npz(itu.FAMILY963_STAGE4)
        itu.assert_size_equals(problem.genotype, 3218, 3)
        itu.assert_problem_stats(problem, 19308, 19286, 23)
        assert_equal(len(list(problem.families(genotyped=False))), 1, 'Incorrect number of families')

        phaser = family_sib_comparison_phaser()
        phaser.run(problem)

        itu.assert_problem_stats(problem, 19308, 19286, 23)

    def test_family_12(self):
        '''Test comparing sibs with non-genotyped parents (stage 4).'''
        problem = io.read_npz(itu.FAMILY12_STAGE2)
        itu.assert_size_equals(problem.genotype, 3218, 7)
        itu.assert_problem_stats(problem, 45052, 42162, 237)
        assert_equal(len(list(problem.families(genotyped=False))), 1, 'Incorrect number of families')

        phaser = family_sib_comparison_phaser()
        phaser.run(problem, PhaseParam(single_member=1))

        itu.assert_problem_stats(problem, 45052, 42162, 237)

    def test_family_2003_need_poo_alignment(self):
        '''Test comparing sibs with non-genotyped parents (stage 4). This case highlights
        the need to align POO-phases, i.e., swap founder haps to correctly patch families at individual
        ID 28412 (our original index 386; in this problem, index 10).'''
        problem = io.read_npz(itu.FAMILY2003_STAGE3)
        itu.assert_size_equals(problem.genotype, 3218, 9)
        itu.assert_problem_stats(problem, 57924, 43339, 85)
        assert_equal(len(list(problem.families(genotyped=False))), 1, 'Incorrect number of families')

        #f = problem.families(genotyped=False)[0]
        #print f.member_list
        #print problem.pedigree.sample_id
        phaser = family_sib_comparison_phaser()
        phaser.run(problem)

        itu.assert_problem_stats(problem, 57924, 57515, 85)

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
