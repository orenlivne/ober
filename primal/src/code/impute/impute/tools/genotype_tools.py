'''
============================================================
Algorithms in related to Pedigree and Genotype objects.

Created on July 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, collections, itertools
from impute.data.constants import MISSING, PATERNAL, MATERNAL, ALLELES, OK, SNP, SAMPLE

####################################################################################
# Genotype Operators, Types and Manipulation
####################################################################################
class GenotypeOperator(object):
    '''Returns a logical array of a logical operation elementwise applied to a list of genotype arrays.'''
    def __init__(self, operator, *g):
        '''operator = operation on (g0f,g0m,...,gnf,gnm).'''
        self.g = g
        self.operator = operator
    
    def __getitem__(self, *args):
        '''Concetenate haplotypes into one argument list (g0f,g0m,...,gnf,gnm) and pass it to operator.'''
        # args = args[0] if len(args) > 1 else args
        slice_args = args[0] if isinstance(args[0], tuple) else (args[0],)
        return self.operator(*(reduce(tuple.__add__, ((g.__getitem__(slice_args + (PATERNAL,)),
                                                       g.__getitem__(slice_args + (MATERNAL,)))
                                                     for g in self.g))))

def action(operator, *g):
    '''A facade that applies an operator to a genotype argument set via slicing. Call as
    apply(operator, g0, ...)[:,:] or [:,x] or [x:y,:], etc.'''
    return GenotypeOperator(operator, *g)

def missing(*g):
    '''Returns a logical array indicating whether the genotype has a missing vlaue.'''
    return np.logical_or(g[PATERNAL] == MISSING, g[MATERNAL] == MISSING)

def not_missing(*g):
    '''Returns a logical array indicating whether the genotype has a missing vlaue.'''
    return np.logical_and(g[PATERNAL] != MISSING, g[MATERNAL] != MISSING)

def homozygous(*g):
    '''Returns a logical array indicating whether the genotype g = (g[0],g[1]) is homozygous.'''
    return np.logical_and(g[PATERNAL] != MISSING, g[MATERNAL] == g[PATERNAL])

def heterozygous(*g):
    '''Returns a logical array indicating whether the genotype g = (g[0],g[1]) is heterozygous.'''
    return np.logical_and(np.logical_and(g[0] != MISSING, g[1] != MISSING), g[0] != g[1])

is_missing = lambda g: action(missing, g)
is_homozygous = lambda g: action(homozygous, g)
is_heterozygous = lambda g: action(heterozygous, g)

def where_missing(g, sample):
    '''Return the SNP indices at which the genotype array g is missing.'''
    return np.where(is_homozygous(g)[:, sample])[0]

def where_homozygous(g, sample):
    '''Return the SNP indices at which the genotype array g is homozygous.'''
    return np.where(is_homozygous(g)[:, sample])[0]

def where_heterozygous(g, sample):
    '''Return the SNP indices at which the genotype array g is heterozygous.'''
    return np.where(is_heterozygous(g)[:, sample])[0]

def index_of(g, value):
    '''Return the index of an allele in a genotype. 'value' is an array of alleles; each allele
    is compared against the corresponding g row. The returned value is the logical array 
    (where(g[0]==value or g[1]==value).'''
    return np.where(np.logical_or(g[:, PATERNAL] == value, g[:, MATERNAL] == value))[0]

def index_first_missing(g):
    '''Return a tuple of (row indices in a genotype array in which the first allele is missing,
    row indices in a genotype array in which the second allele is missing but not the first).'''
    return (np.where(g[:, PATERNAL] == MISSING)[0],
            np.where(np.logical_and(g[:, PATERNAL] != MISSING, g[:, MATERNAL] == MISSING))[0])

def where_heterozygous_and_fraction(h, sample):
    '''Return the location of filled het snps in a haplotype and their fraction out of
    all non-hom-haps.'''  
    het_snps = where_heterozygous(h, sample)
    num_het_snps = len(het_snps)
    return (het_snps, (1.0 * num_het_snps) / (len(where_missing(h, sample)) + num_het_snps)) 

def complete_haplotype(h, g, error_status=OK):
    '''Vectorized version. Given one determined haplotype in h and a genotype g, set the other h entry unless
    the genotype has missing values. So if h=(1,0) and g=(0,2), h will not change, but if
    h=(1,0) and g=(1,2), h will be set to (1,2). Returns the updated haplotype array h.'''
    g_exists = ((g[:, PATERNAL] != MISSING) & (g[:, MATERNAL] != MISSING) & (error_status == OK))
    for p in ALLELES:
        q = 1 - p
        hp_exists = g_exists & (h[:, p] != MISSING)
        fill = hp_exists & (h[:, p] == g[:, p])
        h[fill, q] = g[fill, q]
        fill = hp_exists & (h[:, p] != g[:, p])
        h[fill, q] = g[fill, p]
    return h

def complete_haplotype_partial(h, g):
    '''Same as complete_haplotype(), but does not assume nor check that genotype g exists for all entries.
    g may be partially-called.'''
    h_exists = (h != MISSING)
    
    for p in ALLELES:
        q = 1 - p
        hp, gp, gq = h[:, p], g[:, p], g[:, q]
        hp_exists = h_exists[:, p]
        # Stem cases: (A) h=10, g=11 or 12 (B) h=10, g=20 or 21
        fill = hp_exists & (((hp == gp) & (gp != MISSING)) | ((gp == 0) & (gq == 3 - hp)))  
        h[fill, q] = g[fill, q]
        # Stem case: h=10, g=01 or 02
        fill = hp_exists & (hp != gp) & (gp != 0) & (gp != gq)
        h[fill, q] = g[fill, p]
    return h

# TODO: replace by argument slicing in complete_haplotype()
def complete_haplotype_sample(h, g, snps, sample, error_status=OK):
    '''Vectorized version a single sample. Given one determined haplotype in h and a genotype g, set the other h entry unless
    the genotype has missing values. So if h=(1,0) and g=(0,2), h will not change, but if
    h=(1,0) and g=(1,2), h will be set to (1,2).'''
    g_exists = ((g[:, PATERNAL] != MISSING) & (g[:, MATERNAL] != MISSING) & (error_status == OK))
    for p in ALLELES:
        q = 1 - p
        hp_exists = g_exists & (h[snps, sample, p] != MISSING)
        fill = hp_exists & (h[snps, sample, p] == g[:, p])
        h[snps[fill], sample, q] = g[fill, q]
        fill = hp_exists & (h[snps, sample, p] != g[:, p])
        h[snps[fill], sample, q] = g[fill, p]

# TODO: replace by argument slicing in complete_haplotype()
def complete_haplotype_single(h, g, error_status=OK):
    '''Single-genotype version. Given one determined haplotype in h and a genotype g, set the other h entry unless
    the genotype has missing values. So if h=(1,0) and g=(0,2), h will not change, but if
    h=(1,0) and g=(1,2), h will be set to (1,2).'''
    g_exists = ((g[PATERNAL] != MISSING) & (g[MATERNAL] != MISSING) & (error_status == OK))
    for p in ALLELES:
        q = 1 - p
        hp_exists = g_exists & (h[p] != MISSING)
        if hp_exists and (h[p] == g[p]):
            h[q] = g[q]
        if hp_exists and (h[p] != g[p]):
            h[q] = g[p]

####################################################################################
# Genotype Errors
####################################################################################
def empty_errors_array():
    '''Return an empty errors array.'''
    return np.array([[], []], dtype=int)

def concatenate_errors(a, b):
    '''Concatenate two error arrays.'''
    return np.concatenate((a, b), axis=1)

####################################################################################
# Families
####################################################################################
def genotyped(problem, samples):
    '''Return the set of genotyped members in a collection.'''
    return [x for x in samples if problem.is_genotyped(x)]

def genotyped_members(problem, family):
    '''Return the set of genotyped members in a family.'''
    return [x for x in family.member_list if problem.is_genotyped(x)]

def genotyped_children(problem, family):
    '''Return the set of genotyped children in a family. Could be applied either to a Problem or a Pedigree
    object, since they both support the interface that is called by this method.'''
    return [x for x in family.children_list if problem.is_genotyped(x)]

####################################################################################
# Simple Genotype Statistics, Allele Labeling
####################################################################################
def allele_frequencies_by_snp_old(g):
    '''Return the allele 1,2 frequencies in a genotype array g. Group by SNP.'''
    b = np.where(g == 1)[0]
    c = collections.Counter(b)
    c1 = np.array([c[x] for x in sorted(c.keys())])
    b = np.where(g == 2)[0]
    c = collections.Counter(b)
    c2 = np.array([c[x] for x in sorted(c.keys())])
    c = 1.0 * (c1 + c2)
    return np.array([c1 / c, c2 / c])

def allele_frequencies_by_snp_old2(g):
    '''Return the allele 1,2 frequencies in a genotype array g. Group by SNP. Faster than
    allele_frequencies_by_snp_old() but runs into MemoryError in np.where() for large arrays.'''
    num_snps = g.shape[0]
    count = np.zeros((2, num_snps), dtype=np.uint)
    for a in ALLELES:
        i = np.array([(k, len(list(l))) for k, l in itertools.groupby(np.where(g == a + 1)[0])], dtype=np.uint)
        count[a, i[:, 0]] = i[:, 1]
    return (1.0 * count) / np.sum(count, axis=0)[np.newaxis]

def allele_frequencies_by_snp(g):
    '''Return the allele 1,2 frequencies in a genotype array g. Group by SNP. Does not suffer
    from MemoryErrors.'''
    count = np.array([(len(np.where(g[i] == 1)[0]), len(np.where(g[i] == 2)[0])) for i in xrange(g.shape[0])]).transpose()
    return (1.0 * count) / np.tile(np.sum(count, axis=0)[np.newaxis], (2, 1))

def allele_frequencies(g):
    c1, c2 = (len(np.where(g == 1)[0]), len(np.where(g == 2)[0]))
    c = 1.0 * (c1 + c2)
    return (c1 / c, c2 / c)

def swap_alleles(g):
    '''Swap the allele labels in a genotype array g.'''
    filled = np.where(g != MISSING)
    g[:][filled[SNP], filled[SAMPLE]] = 3 - g[:][filled[SNP], filled[SAMPLE]]

#---------------------------------------------
# Private Methods
#---------------------------------------------
