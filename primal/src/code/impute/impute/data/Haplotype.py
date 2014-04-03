'''
============================================================
A data set of inferred haplotypes of a list of individuals
in a pedigree. This is a 3-D array (individual x SNP x allele)
similar to Haplotype, but here the allele entries are ordered
so that the first is paternally-inherited and the second
is maternally-inherited. 

This is the request passed between phasing stages.
It also holds a reference to the genotype data. Both genotype
and haplotype data arrays may be modified by filters.

We also keep track of quality/confidence scores (QC) in
[0,1] for each haplotype call.

Created on June 13, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, impute as im
from impute.data.Genotype import Genotype

class Haplotype(Genotype):
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, data, snp, sample_id, qc=None, hap_type=None):
        '''Initialize a Haplotype data set from a data array, snp metedata record array corresponding
        to the first dimension elements of data, array of sample IDs corresponding to the second dimension
        elements of data, and optionally:
        - QC measure of each haplotype (float array, entries in [0,1], 0=lowest quality, 1=highest quality)
        - hap_type = array of tags/codes of each haplotype (unphased=0, phased=1, phased with paternal origin=2). Default=1'''
        Genotype.__init__(self, data, snp, sample_id)
        if qc is None:
            # Default: all confidence scores are set to 1 (highest)
            self.qc = np.zeros_like(self.data)
            self.qc[:] = 1.0
        else:
            self.qc = qc

        if hap_type is None:
            # Default: set all haplotype tagging to phased
            self.hap_type = np.zeros(self.data.shape[:2], dtype=np.byte)
            self.hap_type[:] = im.constants.PHASED
        else:
            self.hap_type = hap_type

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return 'Haplotype[snps=%d, samples=%d]' % (self.num_snps, self.num_samples) 
    
    def haplotype(self, person, index):
        '''Return the haplotype of a person for all SNPs (index specifies paternal/maternal/both, if it is a range).'''
        return self.data[person, :, index]

    def paternal_haplotype(self, person):
        '''Return the paternal haplotype of a person.'''
        return self.data[person, :, Haplotype.PATERNAL]

    def maternal_haplotype(self, person):
        '''Return the maternal haplotype of a person.'''
        return self.data[person, :, Haplotype.MATERNAL]
            
    #---------------------------------------------
    # Properties
    #---------------------------------------------
