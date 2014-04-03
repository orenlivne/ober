'''
============================================================
Instantiates data structure types via factory methods.

Created on June 12, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
from impute.data.Genotype import Genotype
from impute.data.Haplotype import Haplotype
from impute.data.Pedigree import Pedigree
from impute.data.GenotypeLazy import GenotypeLazy

####################################################################################
class GenotypeFactory(object):
    '''Instantiates Genotype sub-types.'''
    
    @staticmethod
    def new_instance(clazz, data, snp, sample_id=None, qc=None, hap_type=None, lazy_load=False):
        '''
        Construct a genotype/haplotype set from data arrays:
        - clazz: 'genotype' or 'haplotype' - class type
        - sample_id: genotyped individuals' ID set; if None, using consecutive IDs
        - snp: SNP name array (usually their rs-names)
        - data: a 3-D genotype data array: (individual x SNP x allele)
        
        Optional flags:
        - lazy_load: if True and clazz='genotype', will load into a Genotype impl that lazily
        loads data rows by snp index
        '''
        if sample_id is None: sample_id = np.array(range(Pedigree.START_ID, Pedigree.START_ID + data.shape[1]))
        
        if clazz == 'genotype':
            if lazy_load: return GenotypeLazy(data, snp, sample_id)
            else: return Genotype(data, snp, sample_id)
        elif clazz == 'haplotype': return Haplotype(data, snp, sample_id, qc=qc, hap_type=hap_type)
        else: raise ValueError('Unrecognized Genotype sub-type ' + clazz)

    @staticmethod
    def empty_from_genotype(genotype):
        '''Allocate an empty data array with the same size and metadata as a
        genotype array.'''
        return Haplotype(np.zeros(genotype.shape, dtype=np.byte), genotype.snp, genotype.sample_id)
    
    @staticmethod
    def select(g, clazz, snps):
        '''Return a subset of the SNPs wrapped in a new Haplotype object (a shallow-copy).'''
        gs = GenotypeFactory.new_instance(clazz, g.data[snps, :, :], g.snp[snps], sample_id=g.sample_id)
        gs.map = g.map[snps]
        return gs
    
    @staticmethod
    def select_samples(g, clazz, samples):
        '''Return a subset of the samples wrapped in a new Haplotype object (a shallow-copy).'''
        gs = GenotypeFactory.new_instance(clazz, g.data[:, samples, :], g.snp, sample_id=g.sample_id[samples])
        gs.map = g.map
        return gs
