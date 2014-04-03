'''
============================================================
IBD segment identification using IBS segments + single-locus
IBD posterior probability estimation.

-----------------------------------
OBSOLETE - USE ibd_hmm INSTEAD!
-----------------------------------

Created on December 19, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
from impute.tools import recode

#---------------------------------------------
# Constants
#---------------------------------------------
# T-state (IBS state) corresponding to hash code of (recoded genotype A, recoded genotype B).
# Hash code = 3*(rA-2) + rB-2 = 3*rA + rB - 8 where rA, rB are the recoded genotypes.
# States T=2,6 correspond to IBS=0, the rest to IBS>=1.
__HASH_TO_T = np.array([1, 3, 2, 5, 7, -5, -2, -3, -1])

#---------------------------------------------
# Methods
#---------------------------------------------
####################################################################################
def prob_ibd_single_locus(problem, id1, id2, snps, params):
    '''Is the segment s (defined by the SNP array snps, which may be a frame within the actual segment)
    an IBD segment between samples id1 and id2 or not? Outputs an IBD probability estimate.
    
    This estimate is based on a single-locus IBD posterior probability estimation. Assuming loci
    are statistically independent.'''
    
    # Gather prior information: allele frequencies and condensed identity coefficients
    p = problem.info.allele_frequency(1)
    sample_id = problem.pedigree.sample_id
    _, Delta = params.id_coefs(sample_id[id1], sample_id[id2])
    print 'Delta', Delta
    
    # Load and hash genotype pairs at each SNP
    g = problem.g
    r1, r2 = recode.recode_single_genotype(g[snps, id1, :]), recode.recode_single_genotype(g[snps, id2, :])
    gg_hash = __HASH_TO_T[3 * r1 + r2 - 8]
    
    # Compute posterior (loop over T-states, bulk-set corresponding entries in output array)
    prob_ibd = np.zeros_like(snps, dtype=np.float) 
    for t in __HASH_TO_T:
        index = np.where(gg_hash == t)[0]
        print 'State T=%d, #occurrences %d' % (t, len(index))
        prob_ibd[index] = ibd_posterior_gg(Delta, p[index], t)
    return prob_ibd

def ibd_posterior_hh(f, p):
    '''Haplotype-haplotype IBD posterior probability for a shared allele with frequency p.'''
    return 1 / (1 + p / (f / (1 - f)))

def ibd_posterior_gg(Delta, p, t):
    '''Genotype-Genotype (GG) IBD posterior probability. Returns a 4-array of arrays corresponding
    to the GG T-state t=+-1,+-3,+-5,+-7.'''
    q = 1 - p
    if t < 0:
        # Allele 1 of Genotype A (A1) is 2, reverse alleles' roles 
        p, q, t = q, p, -t
        
    if t == 1:
        alpha = p * (Delta[2] + Delta[4] + Delta[6]) + p ** 2 * Delta[7] + Delta[0]
        beta = p * (Delta[1] + p * (Delta[3] + Delta[5]) + p ** 2 * Delta[8])
    elif t == 2:
        alpha = np.zeros_like(p)
        beta = np.ones_like(p) # Dummy; note that P(IBD|T2)=0
    elif t == 3: 
        alpha = p * Delta[7] + Delta[2]
        beta = 2 * p * (Delta[3] + p * Delta[8])
    elif t == 5:
        alpha = p * Delta[7] + Delta[4]
        beta = 2 * p * (Delta[5] + p * Delta[8])
    elif t == 7:                      
        alpha = np.ones_like(p) * 0.5 * (Delta[6] + 2 * Delta[7])
        beta = 2 * p * q * Delta[8]
    else:
        raise ValueError('Unsupported GG T-state %d' % (t,))
    return alpha / (alpha + beta)

#---------------------------------------------
# Private Methods
#---------------------------------------------
