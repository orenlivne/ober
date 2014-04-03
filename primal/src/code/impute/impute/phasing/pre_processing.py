'''
============================================================
Pre-processing the genotype data: generating SNP annotations
and other useful information stored in ProblemInfo.

In particular, pre-imputation allele frequencies are
computed for use in IBD segment Bayes calculations. They
may be recalculated after imputation, since phasing may
bias the allele frequencies in SNPs with large % of missing
genotypes.

Created on October 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, statutil
from impute.tools import recode
from impute.phasing.phase_core import new_phaser_chain
from chain import Filter

####################################################################################
def __handle_estimate_genotype_frequencies(self, request):
    '''Estimate genotype frequencies from the genotype data and save them in ProblemInfo.'''
    # Load problem fields 
    problem = request.problem
    snp_metadata = problem.info.snp
    snp_count = snp_metadata['count']

    # Recode genotypes to a single number
    r = recode.recode_single_genotype(problem.genotype.data)
    
    # Count genotype appearances for each SNP, and save in SNP annotation array.
    # The frequency table column order matches the GENOTYPE_CODE array. This includes filled
    # and missing genotypes: (1,1),(1,2),(2,2),(0,0).
    for col, genotype_code in enumerate(recode.GENOTYPE_CODE.itervalues()):
        snp_count[:, col] = statutil.hist(np.where(r == genotype_code)[0], problem.num_snps)
    
    # Calculate frequencies
    snp_metadata['frequency'] = statutil.scale_row_sums(snp_count.astype('float'))

    return False

####################################################################################
def __handle_phased_samples(self, request):
    '''Copy data of already-phased samples from the genotype to haplotype array, if the
    PhaseParam object has a non-None selected_samples property. In this scenario, samples
    in selected_samples are to be phased, while the remaining are already phased.'''
    # Load problem fields 
    if request.params.selected_mode:
        phased = np.array(list(set(range(request.problem.num_samples)) - set(request.params.selected_samples)))
        g, h = request.problem.data
        h[:, phased, :] = g[:, phased, :]
    return False

####################################################################################
'''Main pre-processing chain - a facade.'''
def prepare_phaser(next_filter=None, debug=False, print_times=False):
    return new_phaser_chain([
                             estimate_frequencies_phaser,
                             phased_samples_phaser,
                             ], name='Pre-Processing', debug=debug, print_times=print_times)

estimate_frequencies_phaser = Filter(name='* Genotype Freq', handle=__handle_estimate_genotype_frequencies)
phased_samples_phaser = Filter(name='* Phased Samples', handle=__handle_phased_samples)
