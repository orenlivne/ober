'''
Test phasing within nuclear families (misc tests that don't
============================================================
fall under the specific TestPhaseFamily*.py suites).

Created on October 8, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest
from impute import impute_test_util as itu
from numpy.ma.testutils import assert_equal
from impute.phasing.phase_trivial import trivial_phaser
from impute.data import io
from impute.phasing.phase_family import family_child_comparison_phaser, family_phaser
from impute.phasing.phase_core import new_phaser_chain
from impute.tools.param import PhaseParam

class TestPhaseFamily(unittest.TestCase):
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
        self.phaser = new_phaser_chain([trivial_phaser(), family_phaser(), family_child_comparison_phaser()])

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_phase_family(self):
        '''Check phasing trivial cases in all genotyped trios.'''
        problem = io.read_plink(pedigree=itu.HUTT_PED, prefix=itu.GENOTYPE_SAMPLE, haplotype=None)
        itu.assert_size_equals(problem.genotype, 8, 1415)
        assert_equal(len(problem.trios()), 869, 'Unexpected # of genotyped trios')
        self.phaser.run(problem, PhaseParam(debug=False))
        itu.assert_problem_stats(problem, 22640, 20225, 144)

    def test_global_ibd_segments(self):
        '''Test the construction of the global IBD segment dictionary.'''
        problem = io.read_npz(itu.FAMILY225_STAGE1)
        ibd = problem.info.ibd
        assert_equal(ibd.length, 0, 'Unexpected initial # of IBD segments')
        self.phaser.run(problem)
        assert_equal(ibd.length, 21, 'Unexpected final # of IBD segments before grouping')
        ibd.group_to_disjoint()
        assert_equal(ibd.length, 32, 'Unexpected final # of IBD segments after grouping')

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
