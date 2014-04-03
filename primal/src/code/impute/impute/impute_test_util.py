'''
============================================================
Impute code testing utilities.

Created on May 31, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import networkx as nx, numpy as np, test_util as tu
from numpy.testing.utils import assert_equal
from impute.data.problem import Problem
from impute.data import io, io_pedigree, io_genotype
from impute.phasing.examples import OUT_PHASING

#---------------------------------------------
# Constants
#---------------------------------------------
'''Hutterites pedigree - test file name. This file was obtained from Jessica and Gaixin
(minped.3671).'''
HUTT_PED = tu.abs_path('pedigree/hutterites.tfam')
HUTT_GENOTYPED_PED = tu.abs_path('pedigree/hutterites.genotyped.tfam')
'''Small test pedigree'''
SMALL_FILE = tu.abs_path('network/small.tfam')

'''Genotype hutt data with 1 SNP'''
GENOTYPE_SAMPLE = tu.abs_path('sample/sample')

'''Synthetic data on a duo (proband + single parent) with several SNPs'''
GENOTYPE_DUO = tu.abs_path('duo/duo')
# Contains solution to the trio phasing problem: parent genotypes and child haplotypes
GENOTYPE_DUO_SOLUTION = tu.abs_path('duo/duo_solution')

'''Synthetic data on a trio (proband + parents) with several SNPs'''
GENOTYPE_TRIO = tu.abs_path('trio/trio')
# Contains solution to the trio phasing problem: parent genotypes and child haplotypes
GENOTYPE_TRIO_SOLUTION = tu.abs_path('trio/trio_solution')

'''Data on all SNPs and a single nuclear family. If not otherwise specified, data is from chr22.'''
FAMILY7 = tu.abs_path('family7/family7')
FAMILY13 = tu.abs_path('family13/family13')
FAMILY13_STAGE2 = tu.abs_path('family13/family13_stage2.npz')
FAMILY13_STAGE4 = tu.abs_path('family13/family13_stage4.npz')
FAMILY12_STAGE2 = tu.abs_path('family/family12_stage2.npz')
FAMILY225_STAGE1 = tu.abs_path('family/family225_stage1.npz')
FAMILY225_STAGE2 = tu.abs_path('family/family225_stage2.npz')
FAMILY945_ONE_PARENT_STAGE2 = tu.abs_path('family/family945_one_parent_stage2.npz')
FAMILY13_STAGE3_PLINK = tu.abs_path('family13/family13')
FAMILY13_STAGE3 = tu.abs_path('family13/family13_stage3.npz')
# Debugging error that appeared after stage 2 - child had wrong hap color in the middle of long IBD section  
FAMILY4_STAGE1 = tu.abs_path('family/family4_stage1.npz')
FAMILY4_STAGE3 = tu.abs_path('family/family4_stage3.npz')
FAMILY2003_STAGE3 = tu.abs_path('family/family2003_stage3.npz')
# Families of quasi-founder sibs (parents not genotyped) that need POO-alignment
SIB_FOUNDERS_STAGE3 = tu.abs_path('family/sib_founders_stage3.npz')
SIB_FOUNDERS_1049_STAGE3 = tu.abs_path('family/sib_founders_1049.npz')

'''Family (79,324) that Gaixin found to have too many zeroed-out genotypes at SNP 2.'''
FAMILY_TOO_ZEROED_STAGE1 = tu.abs_path('family/family_too_zeroed_stage1.npz')
FAMILY_TOO_ZEROED_STAGE2 = tu.abs_path('family/family_too_zeroed_stage2.npz')
FAMILY_TOO_ZEROED_STAGE4 = tu.abs_path('family/family_too_zeroed_stage4.npz')
FAMILY_TOO_ZEROED_STAGE5 = tu.abs_path('family/family_too_zeroed_stage5.npz')

# Family with one WGS sample (911) that's quite isolated from the rest of the pedigree
CHR18_FAMILY_ISOLATED = tu.abs_path('chr18/family/family_isolated.npz')

'''Distant IBD test sets (node + its pedigree neighbors).'''
NBHRS1298_STAGE4 = tu.abs_path('ibd_distant/nbhrs1298_stage4.npz')
NBHRS1298_STAGE4_WITH_PARENTS = tu.abs_path('ibd_distant/nbhrs1298_stage4_parents.npz')
FAMILY963_STAGE4 = tu.abs_path('family/family963_stage4.npz')

'''GERMLINE IBD test files'''
GERMLINE_HAP = tu.abs_path('ibd_germline/hap.txt')

'''Other IBD files'''
# A limited set of lethal segments - IBDLD IBD-1 segments between parents - chr 22
IBD_LETHAL_COUPLES = tu.abs_path('ibd/couples_chr22.ibd')
SEGMENT_INDEX_CHR22 = tu.abs_path('ibd/index-segments-chr22')

'''Imputation test files'''
IMPUTE_RARE = tu.abs_path('impute/rare.npz')
# A subset of the IMPUTE_RARE sample
IMPUTE_RARE_SAMPLE = tu.abs_path('impute/rare_sample.npz')
# A sample+snp subset of the phased haplotypes in hutt.phased.npz
IMPUTE_PHASED_SAMPLE = tu.abs_path('impute/phased_sample.npz')
# IBD dictionary sample 
IMPUTE_IBD_SAMPLE = tu.abs_path('impute/ibd_sample.txt')

'''IBD Segments in families coloring.'''
# Quasi-founders
QF_FAMILY5_SIBS4 = [250, 597, 631, 1028]
QF_FAMILY5_SIBS4_SEGMENTS = tu.abs_path('color/qf-family5-poo-sibs.txt')
# Non-founder foundation
NF_FAMILY_SIBS = [864, 937, 876, 847]
NF_FAMILY_SEGMENTS = tu.abs_path('color/nf-family.txt')

'''PLINK-imputation-overriding test data directory.'''
OVERRIDE_IMPUTED_BY_PLINK = tu.abs_path('override-imputed-by-plink')

#---------------------------------------------
# External resources that tests rely on.
# TODO: copy to the test resources directory; extract relevant parts to reduce size.
#---------------------------------------------
'''Path to pairwise identity coefficient file.'''
ID_COEF_FILE = OUT_PHASING + '/hutt.id'
'''Path to pairwise kinship file.'''
KINSHIP_FILE = OUT_PHASING + '/hutt.kinship'

#---------------------------------------------
# Methods
#---------------------------------------------
def read_pedigree_from_test_file(file_name, genotyped_id_file=None):
    '''Load a pedigree from a PLINK TFAM file.'''
    data = np.genfromtxt(file_name, np.dtype(int))
    p = io_pedigree.read(file_name, genotyped_id_file=genotyped_id_file)
    assert_equal(p._graph.number_of_nodes(), data.shape[0], 'Incorrect number of nodes')
    assert nx.is_directed_acyclic_graph(p._graph), 'Pedigree is not a DAG'
    return p
    
def assert_segments_almost_equal(segments, expected_data, decimal=6, err_msg='', full_data=False):
    '''Assert that two segment lists are equal up to a tolerance. TODO: replace by Segment classes.'''
    if segments.length != len(expected_data):  # (segments.length > 0) ^ (len(expected_data) > 0):
        raise AssertionError('Wrong segments ' + err_msg + '\nActual Segments:\n' + segments.pprint_segments(True) 
                             + '\nExpected Segments:\n[' + 
                             ',\n '.join('((%-4d, %-4d), (%-8d, %-8d, %7.3f, %d), %s)' % \
                                         (x[0] + x[1] + (repr(x[2]),)) 
                                         for x in expected_data) + ']')    
    try:
        for i in xrange(min(segments.length, len(expected_data))):
            assert_segment_equals_data(segments[i], expected_data[i], full_data=full_data)
    except AssertionError as e:
        raise AssertionError('Segment #' + repr(i) + ': ' + e.message + '\nActual Segments:\n' + segments.pprint_segments(True) 
                             + '\nExpected Segments:\n[' + 
                             ',\n '.join('((%-4d, %-4d), (%-8d, %-8d, %7.3f, %d), %s)' % \
                                         (x[0] + x[1] + (repr(x[2]),)) 
                                         for x in expected_data) + ']')

def assert_size_equals(genotype, num_snps, num_samples):
    '''Assert that a genotype data set's size equals the expected size.'''
    assert_equal(genotype.num_snps, num_snps, 'Incorrect number of SNPS')
    assert_equal(genotype.num_samples, num_samples, 'Incorrect number of samples')

def assert_problem_stats(problem, num_data, num_filled, num_errors, error_rate=0.01):
    '''Assert that a Problem object possesses some desired properties.'''
    h = problem.haplotype
    assert_equal(h.num_data, num_data, 'Incorrect number of total haplotypes')
    assert_equal(h.num_filled, num_filled, 'Incorrect number of filled SNPs')
    # assert_almost_equal(h.fill_fraction, 0.88, 2, 'Incorrect number of initially-filled SNPs')
    assert_equal(problem.num_errors, num_errors, 'Incorrect number of genotype errors')
    assert_equal(problem.num_errors < error_rate * h.num_data, True, '# Genotype errors above expected threshold')

####################################################################################
class Templates:
    '''Hosts sample problems.''' 
    # Sample problem
    PROBLEM_HUT = None
    PROBLEM_FAMILY = None

    @staticmethod 
    def pedigree_hut():
        '''Load the hutterites pedigree.'''
        return read_pedigree_from_test_file(HUTT_PED, genotyped_id_file=GENOTYPE_SAMPLE + '.tfam')
    
    @staticmethod
    def problem_hut():
        '''Load the hutterites data set. Cached since it's large.'''
        if not Templates.PROBLEM_HUT:
            pedigree = Templates.pedigree_hut()
            genotype = io_genotype.read('plink', 'genotype', prefix=GENOTYPE_SAMPLE, load_ids=False)
            Templates.PROBLEM_HUT = Problem(pedigree, genotype) 
        return Templates.PROBLEM_HUT 
        
    @staticmethod
    def haplotype_hut():
        '''Load partially-phased haplotypes for the hutterites sample problem. 
        Cached since it's large.'''
        if not Templates.HAPLOTYPE_HUT:
            Templates.HAPLOTYPE_HUT = io_genotype.read('plink', 'haplotype', tped=GENOTYPE_SAMPLE + '.hap.tped',
                                                       load_ids=False)
        return Templates.HAPLOTYPE_HUT

    @staticmethod
    def problem_family(family, haplotype=False):
        '''Load a single nuclear family data set from PLINK data.
        If haplotype=True, setting initial haplotypes.'''
        return io.read_plink(pedigree=family + '.tfam',
                             prefix=family,
                             haplotype=family + '.hap.tped' if haplotype else None,
                             info=None, idcoef=None, frames=family + '.frm')

#---------------------------------------------
# Private Methods
#---------------------------------------------
def assert_segment_equals_data(segment, data, full_data):
    '''Convert a SegmentSet segment list to a numpy array.''' 
#    for s in segments:
#        print s.snp, s.bp, (s.length, s.num_errors), (hash(Hashable(list(s.samples))),)
    assert_equal(segment.snp, data[0], 'Different SNP range: %s, %s' % (repr(segment.snp), repr(data[0])))
    if full_data:
        assert_equal(segment.bp, data[1][0:2], 'Different base-pair range: %s, %s' % (repr(segment.bp), repr(data[1])))
    assert_equal(sorted(list(segment.samples)), sorted(list(data[2])), 'Different sample set: %s, %s' % (repr(segment.samples), repr(data[2])))
