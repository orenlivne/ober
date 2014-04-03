'''
============================================================
Functions related to recoding genotype pairs as single
numbers.

Created on October 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, itertools as it, util
from impute.data.constants import MISSING, PATERNAL, MATERNAL, INDETERMINATE
from collections import OrderedDict

#---------------------------------------------
# Constants
#---------------------------------------------
# Genotype coding
GENOTYPE_CODE = OrderedDict()
# Include missing genotypes
for genotype in np.concatenate((np.array(list(it.combinations_with_replacement(range(2), 2))) + 1, np.array([[0, 0]])), axis=0):
    GENOTYPE_CODE[tuple(genotype.tolist())] = np.sum(genotype)
# Indices of filled genotypes in the GENOTYPE_CODE array
FILLED_GENOTYPES = np.arange(3)

# A dictionary of (r1,r2)-to-IBS-difference-value, where r1,r2 are recoded values obtained by
# recode_single_genotype(). This is a lookup table,for speed.
IBS_DIFF = {
            (+0, +0): 0,
            (+0, -1): 0,
            (+0, -2): 0,
            (+0, +2): 0,
            (+0, +3): 0,
            (+0, +4): 0,
            
            (-1, +0): 0,
            (-1, -1): 1,
            (-1, -2): 0,
            (-1, +2): 1,
            (-1, +3): 1,
            (-1, +4): 0,
             
            (-2, +0): 0,
            (-2, -1): 0,
            (-2, -2): 1,
            (-2, +2): 0,
            (-2, +3): 1,
            (-2, +4): 1,
             
            (+2, +0): 0,
            (+2, -1): 1,
            (+2, -2): 0,
            (+2, +2): 2,
            (+2, +3): 1,
            (+2, +4): 0,
             
            (+3, +0): 0,
            (+3, -1): 1,
            (+3, -2): 1,
            (+3, +2): 1,
            (+3, +3): 2,
            (+3, +4): 1,

            (+4, +0): 0,
            (+4, -1): 0,
            (+4, -2): 1,
            (+4, +2): 0,
            (+4, +3): 1,
            (+4, +4): 2
            }
_IBS_STATE_FUNC = np.vectorize(lambda x, y: IBS_DIFF[(x, y)])

#---------------------------------------------
# Methods
#---------------------------------------------
def missing_genotype(g):
    '''Return a logical array of indices at which g has at least one missing entry.'''
    return np.min(g, axis=g.ndim - 1) == MISSING

def filled_genotype(g):
    '''Return a logical array of indices at which g has both entries filled.'''
    return np.min(g, axis=g.ndim - 1) != MISSING

def where_missing_genotype(g):
    '''Return the (SNP,sample) indices at which g has at least one missing entry.'''
    return np.where(np.min(g, axis=g.ndim - 1) == MISSING)

def where_filled_genotype(g):
    '''Return the (SNP,sample) indices at which g has both entries filled.'''
    return np.where(np.min(g, axis=g.ndim - 1) != MISSING)

def recode_single_genotype(g):
    '''Recode a genotype data set from (0,0),(0,1),(0,2), (1,1),(1,2),(2,2) to 0,-1,-2, 2,3,4.
    Note that this coding is grouped so that missing data are non-positive. This allows us to compare
    ordered haplotypes with unordered genotypes.'''
    recoded = np.sum(g, axis=g.ndim - 1)
    # Reverse the sign of missing entries from 0,1,2 to 0,-1,-2
    recoded[np.min(g, axis=g.ndim - 1) == MISSING] *= -1 
    return recoded

def where_missing(r):
    '''Return the (SNP,sample) indices where the recoded genotype data array r has both missing entries.'''
    return np.where(r == 0)

def where_called(r):
    '''Return the (SNP,sample) indices where the recoded genotype data array r has both entries
    filled.'''
    return np.where(r > 0)

def where_partial_called(r):
    '''Return the (SNP,sample) indices where the recoded genotype data array r has one filled entry.'''
    return np.where(r < 0)

def where_error(r, r_orig):
    '''Return the (SNP,sample) indices where both the recoded genotype data arrays r, r_orig are called and
    are different.'''
    return np.where((r > 0) & (r_orig > 0) & (r != r_orig))

def where_partial_error(r, r_orig):
    '''Return the (SNP,sample) indices where the recoded genotype data array r has one entry called and
    is incompatible with r_orig.'''
    return np.where(((r == -1) & (r_orig == 4)) | ((r == -2) & (r_orig == 2)))

def where_imputed(r, r_orig):
    '''Return the (SNP,sample) indices where the recoded genotype data array r is either fully or
    partially called and r_orig is not.'''
    return np.where((r != 0) & (r_orig == 0))

def where_full_imputed(r, r_orig):
    '''Return the (SNP,sample) indices where the recoded genotype data array r is fully called and
    r_orig is not fully called.'''
    return np.where((r > 0) & (r_orig <= 0))

def where_partial_imputed(r, r_orig):
    '''Return the (SNP,sample) indices where the recoded genotype data array r is partially called and
    r_orig is not called at all.'''
    return np.where((r < 0) & (r_orig == 0))

def where_still_missing(r, r_orig):
    '''Return the (SNP,sample) indices where the recoded genotype data array r is not called as well as
    r_orig is not fully called.'''
    return np.where((r == 0) & (r_orig <= 0))

def clear_partial_calls(h):
    '''Zero-out partially-called haplotypes.'''
    index = where_partial_called(recode_single_genotype(h))
    h[index[0], index[1], :] = MISSING

def ibs_state(r1, r2):
    '''Return an integer array diff indicating whether the recoded genotype data arrays r1, r2 are
    IBS = 0, 1, or 2.'''
    return _IBS_STATE_FUNC(r1, r2)

def ibs_diff_gh(g, h):
    '''Return a boolean array diff indicating whether the genotype array g and haplotype h array are
    IBS >= 1 (diff=0) or IBS >= 1 (diff=1). Missing values (fully-missing genotype or missing haplotype)
    are treated as indeterminate.'''
    result = ((g[:, PATERNAL] != h) & (g[:, MATERNAL] != h)).astype(np.byte)
    result[(h == MISSING) | ((g[:, PATERNAL] == MISSING) & (g[:, MATERNAL] == MISSING))] = INDETERMINATE
    return result

def concordance_recoded(ra, rb):
    '''Return the concordance between the nx1 singly-recoded genotype arrays a and b. Only takes into
    account entries where both genotypes are fully called.'''
    called = np.where((ra > 0) & (rb > 0))[0]
    ra, rb = ra[called], rb[called]
    concordant = len(np.where(ra == rb)[0])
    discordant = len(np.where(ra != rb)[0])
    return (1.0 * concordant) / (concordant + discordant), concordant, discordant

def concordance_genotypes(a, b):
    '''Return the concordance between the nx2 genotype arrays a and b. Only takes into
    account entries where both genotypes are fully called.'''
    return concordance_recoded(recode_single_genotype(a), recode_single_genotype(b))

#---------------------------------------------
# PLINK- and CGI-related recoding
#---------------------------------------------
# CGI genotype coding
CGI_MISSING_LETTER = 'N'
CGI_GENOTYPES = [x[0] + x[1] for x in list(it.product(CGI_MISSING_LETTER + '01', CGI_MISSING_LETTER + '01'))]

# CGI letter to PLINK recode12 - conversion table
CGI_LETTER_TO_ALLELE = {'N': 0, '0': 1, '1': 2}
ALLELE_TO_CGI_LETTER = util.dict_invert(CGI_LETTER_TO_ALLELE)

'''Convert 0,1,2 allele values in a genotype data array to CGI letters N,0,1.'''
recode_cgi = lambda g: np.vectorize(lambda x: ALLELE_TO_CGI_LETTER[x])(g)

CGI_LETTER_TO_ALLELE_FLIPPED = {'N': 0, '0': 2, '1': 1}
ALLELE_TO_CGI_LETTER_FLIPPED = util.dict_invert(CGI_LETTER_TO_ALLELE_FLIPPED)

'''Same as recode_cgi, only that alleles are also flipped (1->2 and 2->1).'''
recode_cgi_flipped = lambda g: np.vectorize(lambda x: ALLELE_TO_CGI_LETTER_FLIPPED[x])(g)

'''Convert CGI letter genotype to PLINK 1-2 allele coding.'''
recode12 = lambda g: str(CGI_LETTER_TO_ALLELE[g[0]]) + ' ' + str(CGI_LETTER_TO_ALLELE[g[1]])
