'''
============================================================
Centralizes global constants.

Created on May 31, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
#---------------------------------------------
# Constants
#---------------------------------------------
'''Autosomal chromosome numbers.'''
NUM_CHROMOSOMES = 22
CHROMOSOMES = range(1, NUM_CHROMOSOMES + 1)

'''Indicates missing data in adjacency lists, or missing haplotype/genotype record'''
MISSING = 0

'''A flagged genotype error in a Problem error array'''
OK = 0
ERROR = 1

'''Paternally-inherited allele index.'''
PATERNAL = 0
'''Maternally-inherited allele index.'''
MATERNAL = 1
'''Both allele locations'''
ALLELES = (PATERNAL, MATERNAL)
ALLELE_LABELS = {PATERNAL: 'paternal', MATERNAL: 'maternal'}
FLIPPED_ALLELES = (MATERNAL, PATERNAL)
        
'''Trio child index'''
CHILD = 2
'''Haplotype difference involving missing values'''
INDETERMINATE = -1
'''Haplotype is different'''
DIFFERENT = 1

'''Indices into haplotypes and validation Experiment data arrays'''
SNP, SAMPLE = range(2)
'''Haplotype tags'''
UNPHASED, PHASED, PHASED_WITH_ORIGIN, LD_IMPUTED, LD_IMPUTED_WITH_ORIGIN = range(5)

'''No. base pairs in a mega base pair'''
MEGA_BASE_PAIR = 1.0e6

'''Tolerance for float comparison.'''
SMALL_FLOAT = 1e-15
