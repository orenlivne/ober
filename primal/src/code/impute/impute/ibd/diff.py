'''
============================================================
General-purpose Genotypes and haplotype difference measures.

Created on July 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
from impute.data.constants import MISSING, INDETERMINATE
from impute.tools import recode

#---------------------------------------------
# Methods
#---------------------------------------------
def hap_diff(a, b, indeterminate=INDETERMINATE):
    '''For vectors a, b, return the binary difference between two haplotypes. If a or b have missing 
    values, the difference is indeterminate. If a or b is matrix, hap differences are applied
    row-wise.'''
    #if isvector(a) and isvector(b):
    d = np.abs(np.sign(a - b)).astype('int8')
    d[np.where(a == MISSING)] = indeterminate
    d[np.where(b == MISSING)] = indeterminate
    return d

def all_diffs(h, id1, id2, hap1_type=None, hap2_type=None, indeterminate=INDETERMINATE):
    '''Return a list of all haplotype differences between sample indices id1, id2 in the haplotype
    data set haplotype. If hap1_type is specified (PATERNAL/MATERNAL), only one haplotype is
    compared in individual 1 against those of individual 2 (similarly hap2_type).'''
    return np.array([hap_diff(h[:, id1, j], h[:, id2, k], indeterminate=indeterminate) 
                     for j in ((0, 1) if hap1_type is None else (hap1_type,))
                     for k in ((0, 1) if hap2_type is None else (hap2_type,))])

def ibs_state(g, id1, id2):
    '''Return the IBS difference between two haplotypes (IBS=0,1 or 2).'''
    g1, g2 = recode.recode_single_genotype(g[:, id1, :]), recode.recode_single_genotype(g[:, id2, :])
    return recode.ibs_state(g1, g2)

def ibs_diff(g, id1, id2):
    '''Return the IBS difference between two haplotypes (0 if IBS >= 1, 1 if IBS = 0).'''
    g1, g2 = recode.recode_single_genotype(g[:, id1, :]), recode.recode_single_genotype(g[:, id2, :])
    return (recode.ibs_state(g1, g2) == 0).astype(np.byte)

def haps_equal(h, id1, id2, snp_start, snp_stop, hap_comparator):
    '''Return the list of equal haplotypes for a snp_segment=(snp_start,snp_stop)
    in the samples id1,id2.
     
    Haplotype comparison is delegated to self.hap_comparator and must be at least
    self.threshold for a pair to be considered equal.
    
    The function emits tuples (i,j), where i and j are allele indices for which
    hapotypes id1-i and id2-j are equal.'''   
    return [(j, k, hap_comparator(h[snp_start:snp_stop + 1, id1, j], h[snp_start:snp_stop + 1, id2, k])) 
            for j in (0, 1) for k in (0, 1)]

def equal_hap_pairs(h, id1, id2, snp_start, snp_stop, threshold, hap_comparator):
    '''Return the list of equal haplotypes for a snp_segment=(snp_start,snp_stop)
    in the samples id1,id2.
     
    Haplotype comparison is delegated to self.hap_comparator and must be at least
    self.threshold for a pair to be considered equal.
    
    The function emits tuples (i,j), where i and j are allele indices for which
    hapotypes id1-i and id2-j are equal.
    
    hap_comparator(a,b) is a functor that should return a confidence measure in [0,1] of the
    equality of the haplotypes a and b.'''   
    return ((m / 2, np.mod(m, 2)) for m in
             np.where([hap_comparator(h[snp_start:snp_stop + 1, id1, j], h[snp_start:snp_stop + 1, id2, k]) 
                       >= threshold for j in (0, 1) for k in (0, 1)])[0])

def hap_comparator_difference(a, b):
    '''Defines the confidence of equality between haplotypes as the number of equal snps
    out of all snps determined in both.'''
    d = hap_diff(a, b)
    num_filled_snps = np.size(np.where(d != INDETERMINATE))
    return 0.0 if num_filled_snps == 0 else (1.0 * np.size(np.where(d == 0))) / num_filled_snps

def hap_corr_matrix(h, members, hap_type=None):
    '''Return an estimate of the correlation matrix among haps (#equal haps/#determinate haps for each
    pair).'''
    n = len(members)
    # Parametrize all sib pairs
    y, x = np.meshgrid(members, members)
    x, y = x.flatten(), y.flatten()
    print x, y
    if hap_type:
        # All sib pairs and a particular hap[0],hap[1] compbination
        hap_x, hap_y = hap_type
        d = hap_diff(h[:, x, hap_x], h[:, y, hap_y])
        r = np.resize((1.0 * np.sum(d == 0, axis=0)) / np.sum(d != INDETERMINATE, axis=0), (n, n))
    else:
        # All sib pairs and all hap pairs
        d = np.concatenate([hap_diff(h[:, x, j], h[:, y, k]) for j in (0, 1) for k in (0, 1)], axis=1)
        r = np.resize((1.0 * np.sum(d == 0, axis=0)) / np.sum(d != INDETERMINATE, axis=0), (4 * n, 4 * n))
    return r

#---------------------------------------------
# Private Methods
#---------------------------------------------
