'''
============================================================
Similar to Genotype, with lazily-loaded data.

Created on October 7, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import linecache, numpy as np

####################################################################################
class GenotypeLazy(object):
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, data, snp, sample_id):
        '''
        Construct a genotype set from data arrays:
        - snp: SNP metadata record array (contains chromosome, name, morgans, base-pair location)
        - data: name of data file that holds a 3-D genotype data array: (individual x SNP x allele).
        Lazily loaded into a line cache.
        - sample_id: genotyped individuals' ID set
        '''       
        # People's IDs
        self.sample_id = sample_id
        self.num_samples = len(sample_id)
        
        # SNP metadata: SNP label, chromosome number, Genetic distance in Morgans, and
        # base pair location for each SNP
        self.snp = snp
        self.num_snps = len(snp)

        # Lazily load data
        self._data_file = data

    #---------------------------------------------
    # Operators
    #---------------------------------------------
    '''Equality and hashing of objects.'''
    def __key(self): return (self.data.tolist(), self.snp.tolist(), self.sample_id.tolist() if self.sample_id is not None else None)
    def __eq__(self, y): return self.__key() == y.__key()
    def __ne__(self, y): return self.__key() != y.__key()
    def __hash__(self): return hash(self.__key())
    def __repr__(self): return 'GenotypeLazy[snps=%d, samples=%d]' % (self.num_snps, self.num_samples)

    # Methods
    #---------------------------------------------
    #---------------------------------------------
    def data_of_snp(self, index):
        '''Return the genotype data of snp # index. index is 0-based.'''
        # Note: linecache is 1-based
        return np.array(map(np.byte, linecache.getline(self._data_file, index + 1).rstrip('\r\n').rstrip('\n').split(' ')[4:]),
                        dtype=np.byte).reshape([self.num_samples, 2])
