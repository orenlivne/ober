'''
============================================================
Test the end-to-end phasing pipeline and post-processors.

Created on October 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, util, impute as im
from numpy.ma.testutils import assert_equal
from db_gene.snp import mock_dao

class TestPhasePipeline(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        self.problem = im.io.read_plink(pedigree=im.itu.HUTT_PED, prefix=im.itu.GENOTYPE_SAMPLE, haplotype=None)
        im.itu.assert_size_equals(self.problem.genotype, 8, 1415)

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_fill_missing(self):
        '''Testing filling-in missing genotypes.'''
        g = self.problem.genotype
        assert_equal(g.num_filled, 22178, 'Incorrect number of initially-filled genotypes')
        phaser = im.phase_core.new_phaser_chain([im.phase_trivial.trivial_phaser(),
                                                 im.pre_processing.estimate_frequencies_phaser,
                                                 im.post_processing.fill_missing_processor()])
        phaser.run(self.problem)
        im.itu.assert_problem_stats(self.problem, 22640, 18567, 10)
        assert_equal(g.num_filled, 22592, 'Incorrect number of imputed genotypes')

    def test_entire_pipeline(self):
        '''Run the entire pipeline. Using a small # of surrogate parents, for speed.'''
        g = self.problem.genotype
        # Inject a mock DAO so that we don't need the real ID coef file, which is huge here
        self.problem.pedigree._idcoef_dao = mock_dao.IdCoefDao(self.problem.pedigree.num_genotyped)
        phaser = im.phase.build_phasing_pipeline(util.Struct(impute=im.phase.IMPUTE_OPTION.IMPUTE_AND_FILL,
                                                             debug=False, print_times=False, stage=0))
        im.phase.run_phasing_chain(phaser, self.problem, im.PhaseParam(distant_phasing_params=[(0.9, 2, 0.95)]))
        im.itu.assert_problem_stats(self.problem, 22640, 20225, 144)
        assert_equal(g.num_filled, 22640, 'Incorrect number of imputed genotypes')
        assert_equal(g.num_missing, 0, 'Incorrect number of missing genotypes; there should not be any after imputation')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
