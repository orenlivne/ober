'''
============================================================
Post-processing the genotype data: impute using the
phased data; fill in the still-missing entries from a 
multinomial genotype distribution estimated from the data
by the prepare module.  

Allele frequencies are recalculated , since phasing may
have biased the them in SNPs with large % of missing
genotypes.

Created on October 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
from impute.tools import recode
from chain import Filter
from impute.data.constants import SNP, SAMPLE
from statutil import multinomial_elementwise, scale_row_sums
from impute.phasing.phase_core import new_phaser_chain
from impute.tools.recode import FILLED_GENOTYPES
from impute.phasing.pre_processing import estimate_frequencies_phaser
from impute.validation.phasing_validation import Stats

#---------------------------------------------
# Methods
#---------------------------------------------   
def impute_from_fully_called(g, h):
    '''Impute missing genotypes in g from fully-called haplotypes h. The imputation is typically done on the
    genotypes after phasing, which may contain zeroed-out entries found to be Mendelian errors, and thus
    will benefit from imputation here. Partially-filled genotypes are also overridden by haps, if the
    latter are fully called. Returns the number of imputed genotypes'''
    imputed = recode.where_full_imputed(recode.recode_single_genotype(h), recode.recode_single_genotype(g))
    g[imputed[SNP], imputed[SAMPLE], :] = h[imputed[SNP], imputed[SAMPLE], :]
    return len(imputed[0])

#---------------------------------------------
# Request Handlers
#---------------------------------------------   
def __handle_save_stats(self, request):
    '''Save statistics (*before* post-processing) in a request attribute.'''
    request.stats = Stats(request.g_orig, request.problem.haplotype, request.problem.num_errors)
    return False

def __handle_impute_from_fully_called(self, request):
    '''A request handler that imputes missing genotypes in g from fully-called haplotypes h. Clear
    partially-called genotypes.'''
    (g, h) = request.problem.data
    num_imputed = impute_from_fully_called(g, h)
    if request.params.debug: print 'Imputed %d genotypes' % (num_imputed,)
    # Zero-out partial haplotypes. Note: PLINK requires full genotypes, so expect failure
    # if this option is turned off. Therefore it is always on, and we only impute from FULLY called
    # haps above.
    recode.clear_partial_calls(g)
    return False

def __handle_fill_missing_genotypes(self, request):
    '''Fill missing genotype entries by randomly sampling from the multinomial distribution with
    estimated genotype frequencies at the corresponding SNP.'''
    # Load problem fields 
    if request.params.debug:
        print 'Filling missing genotypes from estimated genotype distribution'
    problem = request.problem
    g = problem.genotype.data
    snp_frequency = problem.info.snp['frequency'][:, FILLED_GENOTYPES]

    # Recode genotypes to a single number
    r = recode.recode_single_genotype(g)
    # Find SNP, sample indices of missing data
    missing = recode.where_missing(r)
    
    # Generate random multinomial values; map them to genotype codes
    filled_code = multinomial_elementwise(scale_row_sums(snp_frequency[missing[SNP]])) + 2 
    
    # Fill-in all genotypes of a certain value in a vectorized manner 
    for (genotype, code) in recode.GENOTYPE_CODE.iteritems():
        index = np.where(filled_code == code)[0]
        g[missing[SNP][index], missing[SAMPLE][index], :] = genotype
    return False

####################################################################################
'''Main pre-processing chains - a facade.'''
def save_stats_processor(next_filter=None, debug=False, print_times=False):
    return new_phaser_chain([
                             estimate_frequencies_phaser,
                             Filter(name='* Save stats', handle=__handle_save_stats)
                             ], name='Save Stats', debug=debug, print_times=print_times)

def impute_processor(next_filter=None, debug=False, print_times=False):
    return Filter(name='* Impute from Fully-Called', handle=__handle_impute_from_fully_called)

def fill_missing_processor(next_filter=None, debug=False, print_times=False):
    return Filter(name='* Fill Missing', handle=__handle_fill_missing_genotypes)
