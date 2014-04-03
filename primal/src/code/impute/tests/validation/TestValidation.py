'''
============================================================
Test phasing within nuclear families.

Created on July 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, numpy as np, os, util
from impute.validation import phasing_validation as v
from impute import impute_test_util as itu
from impute.data import io
from numpy.ma.testutils import assert_equal, assert_almost_equal
from impute.phasing.phase_main import main_phaser, nuclear_family_phaser
from impute.plot import plots
from impute.tools import recode

class TestValidation(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        self.problem = io.read_npz(itu.NBHRS1298_STAGE4)
        # Reset seed for deterministic results
        np.random.seed(0)

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_clear_random_portion(self):
        '''Test selecting and clearing a random portion of a genotype set.'''
        # Make sure g has no missing data, easier to check change in g after clearing portion in that case
        g_orig = self.problem.genotype.data.copy()
        g_orig[recode.where_missing_genotype(g_orig)] = 1
        
        g = g_orig.copy()
        g_size = g.shape[0] * g.shape[1]
        assert_equal(len(recode.where_missing_genotype(g)[0]), 0, 'There should be no missing entries in g')
        
        portion = 0.05
        (h, _) = v.clear_random_portion(g, portion)
        
        assert_equal(h.shape[0], int(np.round(portion * g_size)), 'Wrong portion size')
        assert_equal(len(np.where(g != g_orig)[0]), 2 * int(np.round(portion * g_size)),
                     'Wrong change in genotype after clearing portion')

    def test_experiment_statistics(self):
        '''Test methods that return statistics of a genotype set.'''
        fraction = 0.3
        experiment = v.Experiment(self.problem, fraction)
        assert_equal(experiment.total_called, 0, 'Wrong # of initially-called genotypes')
        assert_equal(experiment.total_errors, 0, 'Wrong # of genotype call errors')
        experiment.run(main_phaser())
        assert_almost_equal(experiment.fraction, fraction, 3, 'Wrong fraction of deleted genotypes')
        plots.plot_experiment_stats(experiment)

    def test_parametric_experiment(self):
        '''Test running a validation experiment on a phaser for an array of deleted genotype
        fraction values.'''
        # Make sure g has no missing data, easier to check change in g after clearing portion in that case
        num_experiments = 3
        fraction = np.linspace(0.01, 0.05, num_experiments).astype(np.float32) #+1e-5 # TODO: change to log scale spacing?
        results = v.parametric_experiment(self.problem, fraction, main_phaser(), verbose=False)
        assert_equal(results.shape, [num_experiments], 'Wrong result array size')
        assert_equal(results['deleted_fraction'], fraction, 'Wrong result record')
        plots.plot_parametric_experiment(results, print_stats=False)
        
    def test_phase_and_get_stats(self):
        '''Test gathering phasing statistics in a Stats object.'''
        nuclear_family_phaser().run(self.problem)
        g_orig = self.problem.genotype.data.copy()
        stats = v.Stats(g_orig, self.problem.haplotype, self.problem.num_errors)
        itu.assert_problem_stats(self.problem, 167336, 137693, 102)
        
        # Check stats fields
        partially_filled_genotype = np.where(recode.recode_single_genotype(self.problem.genotype.data) < 0)
        assert_equal(len(partially_filled_genotype[0]), 279, 'Unexpected # of partially-filled genotypes')
        assert_equal(stats.imputed.total, 185, 'Unexpected # of fully-imputed genotypes')
        assert_almost_equal(stats.imputed.fraction, 0.00221, decimal=5, err_msg='Unexpected fraction of fully-imputed genotypes')
        assert_equal(stats.imputed_partial.total, 45, 'Unexpected # of partially-imputed genotypes')
        # Check that printouts don't crash us
        stats.pprint(open(os.devnull, 'wb'))
        # Check that stats is pickable and unpicklable
        out_name = util.temp_filename(suffix='.npz')
        out = open(out_name, 'wb')
        np.savez(out, stats=np.array([stats]))
        out.close()
        loaded = np.load(out_name)
        assert_equal(loaded['stats'][0].imputed_partial.total, stats.imputed_partial.total, 'Pickling & unpickling a Stats object failed')
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
