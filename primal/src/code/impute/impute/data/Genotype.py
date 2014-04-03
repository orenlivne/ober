'''
============================================================
A data set of genotypes of a list of individuals in a
pedigree. This is a 3-D array (individual x SNP x allele).

Created on June 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, util
from util import dict_invert, optimal_insertion_order
from bst import BinarySearchTree
from impute.data import constants

####################################################################################
class Genotype(object):
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, data, snp, sample_id):
        '''
        Construct a genotype set from data arrays:
        - snp: SNP metadata record array (contains chromosome, name, morgans, base-pair location)
        - data: a 3-D genotype data array: (individual x SNP x allele)
        - sample_id: genotyped individuals' ID set
        '''       
        # People's IDs
        self.sample_id = sample_id
        self.data = data
        self._num_snps = self.data.shape[0]
        self._num_samples = self.data.shape[1]
        self._snp_range = None
        
        # SNP metadata: SNP label, chromosome number, Genetic distance in Morgans, and
        # base pair location for each SNP
        self.snp = snp
        # Base-pair-location to snp-index map, lazily-initialized + cached
        base_pair = self.snp['base_pair']
        self._base_pair = base_pair  # np.array([int(base_pair)]) if base_pair.size == 1 else base_pair
        self._bp_to_snp = dict_invert(dict(enumerate(self._base_pair)))
        # Construct a BST for fast bp queries
        self._snp_tree = BinarySearchTree(values=self._base_pair[optimal_insertion_order(self._num_snps)])
        self._snp_index_tree = util.list_index_tree(self._base_pair)
        # A genetic map: lists the two allele letters corresponding to 1 and 2 for each SNP, according
        # their order in the self.snp array.
        self.map = []
        # General metadata, for easy handling of CGI data
        self.metadata = []
                    
        # samples for which the parent-of-origin phase is determined
        self.poo_phase = np.zeros((self._num_samples,), dtype=np.byte)

    @staticmethod
    def empty(snp, sample_id, sz):
        '''Allocate an empty data array of size sz.'''
        return Genotype(snp, sample_id, np.zeros(sz, dtype=np.byte))
    
    def copy(self):
        '''Return a deep copy of this object.'''
        return Genotype(self.data.copy(), self.snp.copy(), self.sample_id.copy())
    
    #---------------------------------------------
    # Operators
    #---------------------------------------------
    def __key(self):
        '''Uniquely-identifying key of this object.'''
        return (self.data.tolist(), self.snp.tolist(),
                self.sample_id.tolist() if self.sample_id is not None else None)

    def __eq__(self, y):
        '''Equality of objects.'''
        return self.__key() == y.__key()

    def __ne__(self, y):
        '''Inequality of objects.'''
        return self.__key() != y.__key()

    def __hash__(self):
        '''Object hash code.'''
        return hash(self.__key())

    def __repr__(self):
        return 'Genotype[snps=%d, samples=%d]' % (self.num_snps, self.num_samples)
     
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def clear(self, snp, person):
        '''Set the data for the individuals whose index is person at snp indices snp
        to MISSING.'''
        self.data[snp, person, :] = constants.MISSING
        
    def nearest_snp(self, bp):
        '''Nearest SNP neighbor to a base pair location bp.'''
        return util.nearest_neighbor_in_list_tree(bp, self._base_pair, self._snp_index_tree)

    def nearest_snp_multiple(self, bp):
        '''Nearest SNP neighbor to each element of a base pair collection bp.'''
        return util.nearest_neighbor_in_list_tree_multiple(bp, self._base_pair, self._snp_index_tree)

    def segment_intersect(self, segment_in_bp):
        '''Return endpoints of the SNP index range contained in the base-pair segment segment.
        If segment does not intersect our SNP range, return None. The range is [start_index,stop_index)
        where start=inclusive and stop is exclsive.'''
        left = self._snp_tree.find_smallest_ge(segment_in_bp[0])
        if left is None: return None
        right = self._snp_tree.find_largest_le(segment_in_bp[1])
        if right is None: return None
        return (self._bp_to_snp[left], self._bp_to_snp[right] + 1)

    def segments_intersect(self, segments_in_bp):
        '''Find SNP range of each segment in the list segments_in_bp whose unit is base pairs.'''
        segments = []
        for x in segments_in_bp:
            y = self.segment_intersect(x)
            if y is not None: segments.append(y)
        return segments
    
    def fill_fraction(self, snps=None, snp_range=None, sample=None, allele=None):
        '''Return the fraction of non-missing data in a sample. If snp_range is specified,
        it takes precedence to the iterable snps. If none is specified, all snps are considered.'''
        
        # TODO: replace by slice iterator
        if allele is not None:
            if sample is not None:
                d = self.data[snp_range[0]:snp_range[1], sample, allele] if snp_range is not None else (self.data[snps, sample, allele] if snps is not None else self.data[:, sample, allele])
            else:
                d = self.data[snp_range[0]:snp_range[1], :     , allele] if snp_range is not None else (self.data[snps, :     , allele] if snps is not None else self.data[:, :, allele])
        else:
            if sample is not None:
                d = self.data[snp_range[0]:snp_range[1], sample, :] if snp_range is not None else (self.data[snps, sample, :] if snps is not None else self.data[:, sample, :])
            else:
                d = self.data[snp_range[0]:snp_range[1], :     , :] if snp_range is not None else (self.data[snps, :     , :] if snps is not None else self.data)
        return 1.0 * np.size(np.where(d != constants.MISSING)[0]) / max(1, np.size(d))
    
    def data_of_snp(self, index):
        '''Return the genotype data of snp # index. index is 0-based.'''
        return self.data[index]

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def shape(self):
        '''Shape of the data array.'''  
        return self.data.shape

    @property
    def num_data(self):
        '''total # data entries.'''  
        return self.data.size
    
    @property
    def num_filled(self):
        '''# filled alleles.'''  
        return len(np.nonzero(self.data)[0])

    @property
    def num_missing(self):
        '''# missing alleles.'''  
        return len(np.where(self.data == constants.MISSING)[0])

# Redundant - cf. fill_fraction()
#    @property
#    def filled_fraction(self):
#        '''Total #filled / # haplotypes.'''  
#        return (1.0 * self.num_filled)/self.num_data

    @property
    def num_snps(self):
        '''#genotyped SNPs for each person.'''  
        return self._num_snps

    @property
    def num_samples(self):
        '''# genotyped people.'''  
        return self._num_samples

    @property
    def base_pair(self):
        '''Base pair locations.'''  
        return self._base_pair

    @property
    def snp_range(self):
        '''Find and cache the range of all SNPs.'''
        if self._snp_range is None: self._snp_range = np.arange(0, self.num_snps)
        return self._snp_range

    @property
    def aligned_samples(self):
        '''Return the list of samples whose parent-of-origin phase has been determined.'''
        return np.where(self.poo_phase)[0]
