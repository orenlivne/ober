'''
============================================================
Reads a segment index from a directory into a class.
The class is presents a facade of uniformly querying across
the entire chromosome, but the underlying implementation
loads only one chromosomal region at one time to save
memory.    

Created on March 19, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, util, itertools as it, os, re
from impute.data.constants import ALLELES, CHROMOSOMES

#---------------------------------------------
# Constants
#---------------------------------------------

#---------------------------------------------
# Methods
#---------------------------------------------

####################################################################################
class SegmentIndex(object):
    '''A segment index. Can be queried at a SNP to find the IBD set that contains
    a certain haplotype.
    
    The Currently-active region is [self._start,self._stop).'''
    
    _REGEX = re.compile('region-(\d+)\.npz')
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, index_dir, chrom=1):
        '''Initialize from a segment index directory.'''
        self._index_dir = index_dir
        self._chrom = 0

        # Force initial load at a test SNP so that the num_samples property doesn't throw an exception
        # For partial index segments that contain only one chromosome during testing,
        # find a chromosome for which a metadata file exists.
        try: 
            found_chrom = it.dropwhile(lambda c: not os.path.exists('%s/chr%d/metadata.npz' % (self._index_dir, c)),
                                       CHROMOSOMES).next()
        except StopIteration:
            raise ValueError('Did not find any metadata under dir %s' % (index_dir,))
        metadata = np.load('%s/chr%d/metadata.npz' % (self._index_dir, found_chrom))
        region_size = metadata['region_size']
        found_region = min(filter(lambda x: x % region_size == 0, map(lambda x: int(x.groups(0)[0]), filter(lambda x: x, map(lambda x: SegmentIndex._REGEX.match(x), os.listdir('%s/chr%d' % (self._index_dir, found_chrom)))))))
        test_snp = found_region
        if chrom > 0 and self._chrom != chrom: self._load_chrom(found_chrom)
        if test_snp < self._start or test_snp >= self._stop: self._load_region(test_snp)
        
        # Detect whether this is a partial index constructed for debugging purposes or a full index
        self.test_index = self._group_index.shape[0] == 1
        # print 'Loaded segment index, found_chrom', found_chrom, 'found_region', found_region, 'test index?', self.test_index

    #---------------------------------------------
    # Operators
    #---------------------------------------------
    def __repr__(self):
        '''Return a pretty-printed-string of a list of segments.'''
        return 'SegmentIndex[index_dir=''%s'']' % (self._index_dir,)

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def num_samples(self):
        '''# samples in index.'''
        return self._group_index.shape[1]
      
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def nearest_left_snp(self, chrom, bp):
        '''Nearest SNP neighbor on the left of a base pair location bp (i.e., largest index snp in the SNP base-pair
        array snp_bp of this object s.t. snp_bp[index] <= bp.'''
        if self._chrom != chrom:
            self._load_chrom(chrom)
        index = util.nearest_neighbor_in_list_tree(bp, self._base_pair, self._snp_index_tree)
        return index if self._base_pair[index] <= bp else index - 1

    def group_index(self, chrom, snp, sample, allele):
        '''Return the IBD group index of a hap at a SNP. Loads a new chromosome and region if necessary,
        so the first call outside the current region will be slower than subsequent calls.'''
        if self._chrom != chrom:
            self._load_chrom(chrom)
        if snp < self._start or snp >= self._stop:
            self._load_region(snp)
    
        if self._region_num >= 0:
            # IBD set number that the haplotype lies in
            relative_snp_index = 0 if self.test_index else snp - self._start
            # Note: IBD group indices are 1-based.
            # print relative_snp_index
            return self._group_index[relative_snp_index, sample, allele]
        else:
            # SNP out of range, so it is not in any IBD group
            return 0 

    def find(self, chrom, snp, sample, allele):
        '''Return an array of haplotypes that are IBD with the haplotype (sample,allele) at SNP index snp,
        including the haplotype itself. If the haplotype is not IBD with anyone else, returns an empty array.'''
        # print relative_snp_index, sample, allele, group
        # print 'find(', chrom, snp, sample, allele, ')'
        group_index = self.group_index(chrom, snp, sample, allele)
        # print snp, self._start, group_index, len(self._groups[snp - self._start])
        group = self._groups[0 if self.test_index else snp - self._start][group_index] if group_index else np.zeros((0, 2), dtype=np.uint)
        return group
        # if group.ndim == 1: group = group[np.newaxis]
        # return group

    def covered_by_training(self, chrom, snp, T):
        '''Return set % of haps that can be imputed from a training sample set T.'''
        T_haps = set(list(it.product(T, ALLELES)))
        return set(tuple(x) for hap in T_haps for x in self.find(chrom, snp, hap[0], hap[1]))

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def _load_chrom(self, chrom):
        '''Load index of a new chromosome number, chrom.'''
        # Load index metadata
        self._chrom = chrom
        metadata = np.load('%s/chr%d/metadata.npz' % (self._index_dir, chrom))
        self._snp = metadata['snp']
        self._region_size = metadata['region_size']
        # Currently-active region is [start,stop)
        self._start = 0
        self._stop = 0 
        self._region_num = -1
        
        # Construct a BST for fast queries of the left-neighboring SNP of a base-pair position
        base_pair = self._snp['base_pair']
        self._base_pair = np.array([int(base_pair)]) if base_pair.size == 1 else base_pair
        self._snp_index_tree = util.list_index_tree(self._base_pair)

    def _load_region(self, snp):
        '''Load regional index around SNP index ''snp''.'''
        r, num_snps = self._region_size, len(self._snp)
        if snp < 0:
            # SNP out of range, to the left of first region.
            # Set convenient fictitious values for active region so that we don't call this
            # function multiple times when find()-ing multiple out-of-range SNPs.
            self._start, self._stop, self._region_num = -1, 0, -1
        elif snp >= num_snps:
            # SNP out of range, to the right of last region.
            # Set convenient fictitious values for active region so that we don't call this
            # function multiple times when find()-ing multiple out-of-range SNPs.
            self._start, self._stop, self._region_num = num_snps, num_snps + 1, -1
        else:
            # SNP in range covered by a region. Load it.
            region_num = snp / r
            self._start, self._stop, self._region_num = region_num * r, min((region_num + 1) * r, num_snps), region_num
            a = np.load('%s/chr%d/region-%d.npz' % (self._index_dir, self._chrom, self._start))
            self._groups = a['groups']
            self._group_index = a['group_index']
    
    def _load_region_by_number(self, region_num):
        '''Load regional index of region number region_num.'''
        r, num_snps = self._region_size, len(self._snp)
        print 'Loading region %d' % (region_num,)
        self._start, self._stop, self._region_num = region_num * r, min((region_num + 1) * r, num_snps), region_num
        a = np.load('%s/chr%d/region-%d.npz' % (self._index_dir, self._chrom, self._start))
        self._groups = a['groups']
        self._group_index = a['group_index']

#    def training_coverage(self, snp, T):
#        '''Return the % of haps that can be imputed from a training sample set T.'''
#        imputable_haps = self.covered_by_training(snp, T)
#        call_rate_hap = (100.*len(imputable_haps)) / num_haps
#        call_rate_genotype = (100.*sum(1 for sample in xrange(info.num_samples) 
#                                       if (sample, PATERNAL) in imputable_haps and
#                                       (sample, MATERNAL) in imputable_haps)) / self._info.num_samples
#        return (call_rate_hap, call_rate_genotype)
