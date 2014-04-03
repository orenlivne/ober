'''
============================================================
Test phasing algorithm for trivial Mendelian cases in all
trios and all SNPs.

Created on July 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest
from impute import impute_test_util as itu
from numpy.ma.testutils import assert_equal
from impute.phasing.phase_trivial import trivial_phaser
from impute.data import io
from impute.phasing.phase_core import new_phaser_chain

class TestPhaseTrivialAll(unittest.TestCase):
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
        self.problem = io.read_plink(pedigree=itu.HUTT_PED, prefix=itu.GENOTYPE_SAMPLE, haplotype=None)
        self.phaser = new_phaser_chain([trivial_phaser()])

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_phase_trivial_cases_all_trios(self):
        '''Check phasing trivial cases in all genotyped trios.'''
        itu.assert_size_equals(self.problem.genotype, 8, 1415)
        assert_equal(len(self.problem.trios()), 869, 'Unexpected # of genotyped trios')
        self.phaser.run(self.problem)
        itu.assert_problem_stats(self.problem, 22640, 18567, 10)
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------