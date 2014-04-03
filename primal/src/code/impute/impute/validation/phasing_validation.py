'''
============================================================
Phasing validation functions: delete a random portion of
the original data set, phase the data, and compare the
generated (imputed) haplotypes with the original data. 

Created on August 30, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, time, util, statutil, numpy as np, collections  # @UnresolvedImport
from impute.data.constants import MISSING
from impute.data.problem import Problem
from impute.tools import recode
from impute.data.constants import SNP, SAMPLE

#---------------------------------------------
# Constants
#---------------------------------------------

#---------------------------------------------
# Methods
#---------------------------------------------
def random_portion(g, n):
    '''Return a random sub-array of size n of a (SNPs,samples) in a genotype data set g. Only filled
    genotypes are considered for inclusion in the sample.'''
    i = recode.where_filled_genotype(g)
    j = statutil.random_index(len(i[0]), n)
    return (i[0][j], i[1][j])

def clear_random_portion(g, fraction):
    '''Set a random sub-array of size fraction*(original size) of a genotype data set g to MISSING, and
    return the original values of the sub-array. If a second output argument is sought, it returns the index
    vector of the sub-array into g.'''
    return clear_index(g, random_portion(g, int(round(fraction * g.shape[0] * g.shape[1]))))

def clear_index(g, index):
    '''Set a random sub-array index of a genotype data set g to MISSING, and return the original
    values of the sub-array. If a second output argument is sought, it returns the index vector
    of the sub-array into g.'''
    portion = g[index[0], index[1], :]
    g[index[0], index[1], :] = MISSING
    return portion, index

def parametric_experiment(problem, deleted_fraction, phaser, params=None, verbose=True):
    '''A parametric validation experiment: run an phaser Experiment for an array of deleted genotype 
    fraction values. Return a record array holding run statistics.'''
    pe = ParametricExperiment(problem)
    return pe.run(deleted_fraction, phaser, params=None, verbose=verbose)

####################################################################################
class StatsBatch(object):
    '''Phasing validation statistics accumulated from a batch run (wraps a 2-D array of chromosome and
    chromosome part Stats objects).'''

    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, stats_array):
        '''Initialize phasing statistics object for the original genotype g/and haplotype set problem.haplotype.'''
        self.stats_array = stats_array
        
    @staticmethod
    def new_instance(num_chroms):
        '''A constructor of an empty stats object to hold num_chromosomes chromosomes.'''
        stats_array = np.empty((num_chroms, 0), dtype=Stats)
        return StatsBatch(stats_array)
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return 'StatsArray[%s]' % (repr(self.stats_array)) 
    
    def add_part(self, chrom, stat):
        '''Add a chromosome part Stats object to the stats data array of chromosome # chromosome.'''
        raise ValueError('to be implemented')
    
    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def stats(self):
        '''Return a tuple containing all experiment statistics: fraction, all call rates, run time.'''
        return (self.called.fraction, self.partial_called.fraction,
                self.errors.fraction, self.errors_partial.fraction)
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------

####################################################################################
class Stats(object):
    '''Phasing validation statistics.'''

    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, g_orig, haplotype, num_errors):
        '''Initialize phasing statistics object for the original genotype g/and haplotype set problem.haplotype.'''
        h = haplotype.data
        r_orig = recode.recode_single_genotype(g_orig)
        r = recode.recode_single_genotype(h)
        
        # Sizes
        self.time = 0
        self.num_snps = haplotype.num_snps
        self.num_samples = haplotype.num_samples
        self.num_genotypes = haplotype.num_data / 2
        self.num_haplotypes = haplotype.num_data
        
        # Arrays
        self.fill = np.array([haplotype.fill_fraction(sample=x) for x in xrange(haplotype.num_samples)])
        
        # Fields
        # A Field factory method
        field = lambda index: StatsField(self, h, index)
        self.called_orig = field(recode.where_called(r_orig))
        self.imputed = field(recode.where_full_imputed(r, r_orig))
        self.imputed_partial = field(recode.where_partial_imputed(r, r_orig))
        self.errors = field(recode.where_error(r, r_orig))
        self.errors_partial = field(recode.where_partial_error(r, r_orig))
        self.called = field(recode.where_called(r))
        self.partial_called = field(recode.where_partial_called(r))
        self.still_missing = field(recode.where_still_missing(r, r_orig))

        # Scalars
        self.num_filled_haplotypes = haplotype.num_filled
        self.num_errors = num_errors  # Redundant
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return 'Stats[%s]' % (repr(self.stats)) 
    
    def pprint(self, out=sys.stdout):
        '''Pretty-Print phasing summary statistics of the updated problem object.'''
        parts = collections.OrderedDict([
                 ('filled', ('Filled haps', self.num_filled_haplotypes)),
                 ('errors', ('Genotype errors', self.num_errors)),
                 ('unphased', ('Remaining unphased', self.num_haplotypes - self.num_errors - self.num_filled_haplotypes)),
                 ('imputed', ('Imputed genotypes', self.imputed.total)),
                 ])        
        out.write('%-20s %-9d (%d snps x %d samples)\n' % \
        ('Total haplotypes', self.num_haplotypes, self.num_snps, self.num_samples,))
        for (label, num) in parts.itervalues():
            out.write('%-20s %-9d %7.3f%%\n' % (label, num, 100.*num / self.num_haplotypes))
        out.write('%-20s %-9.1f (%.1e sec / data)\n' % ('Time [sec]', self.time, 1.*self.time / self.num_haplotypes,))
    
    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def stats(self):
        '''Return a tuple containing all experiment statistics: fraction, all call rates, run time.'''
        return (self.called.fraction, self.partial_called.fraction,
                self.errors.fraction, self.errors_partial.fraction)
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    
####################################################################################
class StatsField:
    '''Holds statistics on a single field/aspect of phasing.'''
    def __init__(self, stats, g, index):
        '''Create a struct of standard properties for one field (#imputed, #genotypes, etc.)
        for the (sample, snp) index array index.'''
        # Owning object
        self.__num_genotypes = stats.num_genotypes
        self.total = len(index[0])
        # Dictionary that represents cut by sample
        self.by_sample = StatsField.__group_by_field(g, index, SAMPLE)

    @property
    def fraction(self):
        '''Return the total % of this property.''' 
        return (1.0 * self.total) / self.__num_genotypes

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    @staticmethod
    def __group_by_field(g, i, field):
        '''Group a test index subset i by field (SNP=0, sample=1).'''
        group_count = util.dict_to_array(statutil.group_by_value(i[field]))
        result = np.zeros((g.shape[field],), dtype=int)
        result[group_count['k']] = group_count['v']
        return result

####################################################################################
class Experiment(object):
    '''A validation experiment: start with a Problem object, clear a certain portion of the
    data, run a phaser, and cross-check the hap results against the original genotype data.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, problem, fraction=None, test_index=None):
        '''Initialize an experiment to be run on a problem, clearing out 'fraction' of the data. If test_index
        is specified, these specific test indices are used; otherwise a random fraction is generated.
        
        If test_index = 'hap', data is read from problem.h (haplotype array). The entire array
        is considered as a test array, but nothing is zeroed out. Useful for phasing result stats.'''
        # Create a working copy of the problem. Only the data is copied.
        if not (fraction is not None) ^ (test_index is not None):
            raise ValueError('Must specify fraction or test_index')
        self.problem = Problem(problem.pedigree, problem.genotype.copy())
        self.h = self.problem.h
        
        # Create test set; save original genotypes in g_orig
        if test_index is None:
            self.fraction = fraction
            self.g_orig, i = clear_random_portion(self.problem.genotype.data, fraction)
        elif test_index == 'hap':
            # Don't clear anything; call everything a test index.
            h = problem.h
            i = tuple(util.flattened_meshgrid(range(h.shape[0]), range(h.shape[1])))
            self.g_orig = problem.g
            self.h = h
            self.fraction = 1.0
        else:
            self.g_orig, i = clear_index(self.problem.g, test_index)
            self.fraction = (1.0 * i[0].size) / (self.h.shape[0] * self.h.shape[1])
        self.num_tests = i[0].size
        self.test_index = i
        self.r_orig = recode.recode_single_genotype(self.g_orig)
        self.fill = self.problem.fill_fraction()[:, SAMPLE]
        self.__recode_single_genotype = None
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return 'Experiment[%s, fraction=%.2f%%]' % (repr(self.problem), self.fraction) 
    
    def run(self, phaser, params=None):
        '''Run phaser (or more generally, a handler) on a problem.'''
        phaser.run(self.problem, params=params)
        self.fill = self.problem.fill_fraction()[:, 1]
    
    def num_test_genotypes(self, field):
        '''Return the number of genotypes in which both alleles were called, broken by field (SNP=0, sample=1).'''
        return self.__group_by_field(np.arange(len(self.test_index[0])), field)
    
    def where_called(self):
        '''Return the indices of genotypes in which both alleles were called.''' 
        # Positive entries of r = called entries
        return recode.where_called(self.recoded_genotype)[0]

    def called(self, field):
        '''Return the number of genotypes in which both alleles were called, broken by field (SNP=0, sample=1).'''
        return self.__group_by_field(self.where_called(), field)

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def test_orig(self):
        '''Return the original set of deleted test genotypes.'''
        i = self.test_index
        return recode.recode_single_genotype(self.h[i[SNP], i[SAMPLE], :])

    @property
    def test_called(self):
        '''Return the called set of haplotypes corersponding to the test genotypes.'''
        i = self.test_index
        return self.h[i[SNP], i[SAMPLE], :]
    
    @property
    def recoded_genotype(self):
        '''Return the genotype test set, recoded as a single number of allele pair.'''
        if self.__recode_single_genotype is None:
            self.__recode_single_genotype = recode.recode_single_genotype(self.test_called) 
        return self.__recode_single_genotype

    @property    
    def total_called(self):
        '''Return the number of genotypes in which both alleles were called.''' 
        return len(self.where_called())

    @property    
    def total_partial_called(self):
        '''Return the number of genotypes in which one alleles was called.''' 
        # Positive entries of r = called entries
        return len(recode.where_partial_called(self.recoded_genotype)[0])
    
    @property    
    def total_errors(self):
        '''Return the number of genotypes that were called incorrectly. (A genotype is an allele pair.)''' 
        # Count entries that were called (positive AND are different than the corresponding original value
        return len(recode.where_error(self.recoded_genotype, self.r_orig)[0])
    
    @property    
    def total_partial_errors(self):
        '''Return the number of genotypes that were called incorrectly. (A genotype is an allele pair.)''' 
        # Count entries that were called (positive AND are different than the corresponding original value
        # This happens when hap=(0,1), genotype=(2,2) or hap(0,2), genotype=(1,1)
        return len(recode.where_partial_error(self.recoded_genotype, self.r_orig)[0])

    @property
    def full_call_fraction(self):
        '''Return the % of correctly fully-called test genotypes.''' 
        return (1.0 * self.total_called) / self.num_tests
    
    @property
    def partial_call_fraction(self):
        '''Return the % of erroneously half-called test genotypes.''' 
        return (1.0 * self.total_partial_called) / self.num_tests

    @property
    def full_error_fraction(self):
        '''Return the % of erroneously fully-called test genotypes.''' 
        return (1.0 * self.total_errors) / self.num_tests

    @property
    def partial_error_fraction(self):
        '''Return the % of erroneously half-called test genotypes.''' 
        return (1.0 * self.total_partial_errors) / self.num_tests

    @property
    def stats(self):
        '''Return a tuple containing all experiment statistics: fraction, all call rates, run time.'''
        return (self.fraction,
                self.full_call_fraction,
                self.partial_call_fraction,
                self.full_error_fraction,
                self.partial_error_fraction)
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __group_by_field(self, i, field):
        '''Group a test index subset i by field (SNP=0, sample=1).'''
        size = self.problem.genotype.data.shape[field]
        group_count = util.dict_to_array(statutil.group_by_value(self.test_index[field][i]))
        result = np.zeros((size,), dtype=int)
        result[group_count['k']] = group_count['v']
        return result

####################################################################################
class ParametricExperiment:
    '''A helper class -- parametric validation experiment.'''
    
    def __init__(self, problem):
        '''Initialize a parametric experiment for a data set.'''
        self.problem = problem
        
    def run(self, deleted_fraction, phaser, params, verbose=False):
        '''Run a parametric experiment for each deleted genotyped % value in deleted_fraction.'''
        stats = np.zeros((len(deleted_fraction),),
                         dtype=[('deleted_fraction', 'f4'),
                                ('full_call_fraction', 'f4'),
                                ('partial_call_fraction', 'f4'),
                                ('full_error_fraction', 'f4'),
                                ('partial_error_fraction', 'f4'),
                                ('run_time_sec', 'f4')])
        a = [x for x in self.__run_experiments(deleted_fraction, phaser, params, verbose)]
        stats[:] = a
        return stats

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __run_experiments(self, deleted_fraction, phaser, params, verbose):
        '''A generator of experiment results.'''
        for fraction in deleted_fraction:
            if verbose:
                print 'Running experiment with %.2f%% deleted genotypes' % (100 * fraction,)
            experiment = Experiment(self.problem, fraction)
            start = time.time()
            experiment.run(phaser, params)
            time_elapsed = time.time() - start
            stats = experiment.stats + (time_elapsed,)
            if verbose:
                print 'Deleted %.2f%%   Call Rate %.1f (%.1f)   Error Rate %.4f (%.4f)   Time [s] %.1f' % \
                (100 * stats[0], 100 * stats[1], 100 * stats[2], 100 * stats[3], 100 * stats[4], stats[5])
            yield stats
