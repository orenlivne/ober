#!/usr/bin/env python
'''
============================================================
Imputation data structures. Read and maintain the imputation
training set (Hutt 98) data and imputation output.
     
Created on December 3, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, impute as im

####################################################################################
class ImputationSet(object):
    '''Reads and maintains the imputation training set (Hutt 98) data and imputation
    output.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, pedigree, genotype, imputed=None, sample_index=None):
        '''Read training set data from in_file. Use the pedigree ''pedigree'' to translate IDs
        between the training set and the phasing set.'''
        # One underscore = protected member
        self.genotype = genotype
        self.pedigree = pedigree
        self.__sample_index = sample_index if sample_index is not None else np.array([pedigree.node_of[x] for x in self.genotype.sample_id])
        self.sample_index_of_typed = dict(zip(self.__sample_index, xrange(genotype.data.shape[0])))
        self.__imputed = imputed if imputed else \
        im.factory.GenotypeFactory.new_instance('haplotype',
                                                np.zeros((self.genotype.num_snps, self.pedigree.num_genotyped, 2), dtype=np.byte),
                                                self.genotype.snp, self.genotype.sample_id,
                                                qc=im.constants.MISSING)

    @staticmethod
    def from_file(pedigree, in_file, imputed=None):
        '''Read training set data from in_file. Use the pedigree ''pedigree'' to translate IDs
        between the training set and the phasing set.'''
        return ImputationSet(pedigree, im.io_genotype.read('npz', 'genotype', file=in_file), imputed=imputed)

    @staticmethod
    def from_problem(problem):
        '''Read a Problem''s genotype data into the genotype array of this object. Read the
        Problem''s haplotype data into the imputed array of this object.'''
        return ImputationSet(problem.pedigree, problem.genotype,
                             sample_index=np.empty((0,)),
                             imputed=problem.haplotype)
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    @staticmethod
    def load(npz_file):
        '''Deserialize to object from an npz file.'''
        data = np.load(npz_file)
        t = im.factory.GenotypeFactory.new_instance('genotype',
                                                    data['t_data'], data['t_snp'], data['t_sample_id'])
        t.map = data['t_map']
        imputed = im.factory.GenotypeFactory.new_instance('haplotype', data['imputed_data'], t.snp, t.sample_id,
                                                          hap_type=data['imputed_hap_type'] if 'imputed_hap_type' in data.files else None)
        return ImputationSet(data['pedigree'][0], genotype=t, imputed=imputed)
        
    def save(self, npz_file):
        '''Serialize to object to an npz file.'''
        np.savez(npz_file,
                 sample_index=self.sample_index,
                 pedigree=np.array([self.pedigree]),
                 imputed_data=self.imputed_data,
                 imputed_hap_type=self.imputed.hap_type,
                 t_data=self.genotype.data,
                 t_snp=self.genotype.snp,
                 t_sample_id=self.genotype.sample_id,
                 t_map=self.genotype.map)

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def num_samples(self):
        '''Return the total number of samples: to-impute + training (# samples that are affy-genotyped).'''
        return self.pedigree.num_genotyped
    
    @property
    def sample_index_to_impute(self):
        '''Return the list of samples to impute, which are not in sample_index.'''
        return np.array(list(set(range(self.pedigree.num_genotyped)) - set(self.sample_index)))

    @property
    def imputed(self):
        '''Return the imputed haplotype object.'''
        return self.__imputed
    
    @property
    def imputed_data(self):
        '''Return the imputed haplotypes data array.'''
        return self.__imputed.data
    
    @imputed_data.setter
    def imputed_data(self, data):
        '''Set the imputed haplotypes data array to the ''data'' array.'''
        self.__imputed.data = data
    
    @property
    def imputed_hap_type(self):
        '''Return the imputed haplotypes tag (phased/unphased) data array.'''
        return self.__imputed.hap_type
    
    @imputed_hap_type.setter
    def imputed_hap_type(self, hap_type):
        '''Set the imputed haplotypes tag (phased/unphased) data array to the ''hap_type'' array.'''
        self.__imputed.hap_type = hap_type

    @property
    def training_data(self):
        '''Return the training genotypes.'''
        return self.genotype.data

    @property
    def sample_index(self):
        '''Return the array of phasing IDs of the training samples.'''
        return self.__sample_index

    @property
    def sample_index_set(self):
        '''The same as training_sample(), but returns a set object.'''
        return set(self.__sample_index)
    
    @property
    def num_snps(self):
        '''Number of SNPs to impute.''' 
        return self.genotype.num_snps
    
    @property
    def snp(self):
        '''SNPs-to-impute metadata array.''' 
        return self.genotype.snp
