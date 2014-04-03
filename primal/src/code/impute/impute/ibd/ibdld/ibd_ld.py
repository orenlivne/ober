'''
============================================================
IBDLD segment retrieval and confidence measures.

Created on July 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import csv, linecache, impute as im
import numpy as np
from impute.ibd import diff
from impute.data.constants import INDETERMINATE
from impute.tools.param import PhaseParam

#---------------------------------------------
# Constants
#---------------------------------------------
'''Lide Han uses this value as a posterior HMM confidence probability cutoff.'''
IBDLD_THRESHOLD = 0.9 

#---------------------------------------------
# Methods
#---------------------------------------------
def ibd_ld_confidence(x, M):
    '''Posterior probability that x is in an IBD segment. Here x is SNP offset from
    the middle of interval [-0.5*M,0.5*M] and M is the length of IBD segment. In practice,
    log(ibd_ld) should be used for numerical stability.'''
    xc = __xcritical(M)
    return IBDLD_THRESHOLD + (1 - IBDLD_THRESHOLD) * (1 - np.exp(-abs((x - 0.5 * M) / xc))) * (1 - np.exp(-abs((x + 0.5 * M) / xc)))

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __xcritical(M):
    '''Critical value of confidence dropping. Large enough for M ~ 1, linear in M for large M.'''
    return 1 + 0.0625 * M / 10

####################################################################################
class IbdSegmentComputerFactory(object):
    '''An abstract factory of an IBD computation scheme - using Mark & Lide's IBDLD.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, file_name):
        '''Initialize an IBDLD cache from an IBDLD output file.'''
        self._cache = IbdSegmentGlobalCacheIbdld(file_name)

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def segment_computer(self, haplotype, **kwargs):
        '''Return the IBD segment cache.'''
        sample_id = kwargs.get('sample_id', None)
        chrom = kwargs.get('chrom', None)
        samples = kwargs.get('samples', None)
        threshold = kwargs.get('threshold', IBDLD_THRESHOLD)
        if any(x is None for x in [chrom, samples]):
            raise ValueError('Must specify sample_id, chrom (chromosome #) and samples')
        return IbdSegmentComputerIbdld(self._cache, haplotype, chrom, samples, sample_id, threshold) 

####################################################################################
class IbdSegmentComputerIbdld(object):
    '''Maintains a cache of IBDLD results of IBD-1 sharing between genotyped samples
    for a single chrom.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, cache, haplotype, chrom, samples, sample_id, threshold, params=PhaseParam()):
        self.__cache = cache 
        self.__sample_id = sample_id
        self.__haplotype = haplotype
        self.__chrom = chrom
        self.__samples = samples
        self.__num_snps = haplotype.num_snps
        self.__threshold = threshold 
        self.__hap_comparator = _hap_comparator_ibdld
        self.__data = haplotype.data
        self.__params = params
        self._bp = haplotype.snp['base_pair']
        self.__num_snps = haplotype.num_snps
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def segments(self):
        '''Return the list of IBD segments shared by individuals id1, id2.'''
        for a in self.__samples:
            a_id = self.__sample_id[a]
            for b in self.__samples:
                b_id = self.__sample_id[b]
                # Prevent duplicate pairs
                if a_id < b_id:
                    segments_in_bp = self.__cache.segments(self.__chrom, a_id, b_id)
                    if self.__params.debug:
                        print segments_in_bp
                    for segment in self.__haplotype.segments_intersect(segments_in_bp):
                        # TODO: find the equal haplotype pairs instead of emitting INDETERMINATE.
                        # It can be DETRIMENTAL and MISLEADING to  try to determine equal haplotypes
                        # from genotypes if only a short interval of non-missing difference snps is
                        # available. See 9-AUG-12 note in Oren's Hutterites binder for an example.
                        pairs = self.__equal_hap_pairs(a, b, segment[0], segment[1])
                        (snp_start, snp_stop) = segment
                        
                        # Debugging printouts                        
                        if self.__params.debug:
                            pairs = list(pairs)
                            print 'Sample indices', a, b, 'segment', segment
                            for j in (0, 1):
                                for k in (0, 1):
                                    ha = self.__data[snp_start:snp_stop + 1, a, j]
                                    hb = self.__data[snp_start:snp_stop + 1, b, k]
                                    confidence = self.__hap_comparator(ha, hb)
                                    print '\t\tpair', (j, k), 'confidence', confidence 
                            for p in pairs:
                                print '\t\tEqual pairs', p

                        for p in pairs:
                            # Note 1: setting bp of segment to last SNP's bp if segment extends beyond
                            # the last SNP.
                            # Note 2: if segment is too short for the SNP grid resolution (is entirely
                            # contained in the gap between consecutive SNPs), ignore it
                            if self.__params.debug:
                                print (a, p[0]), (b, p[1]), snp_start, snp_stop
                            if snp_stop > snp_start:
                                yield im.segment.Segment(segment, ((a, p[0]), (b, p[1])),
                                                         (self._bp[snp_start], im.segment.stop_bp(self._bp, snp_stop, self.__num_snps)))
                            
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __equal_hap_pairs(self, id1, id2, snp_start, snp_stop):
        '''Return the list of equal haplotypes for a snp_segment=(snp_start,snp_stop)
        in the samples id1,id2.
         
        Haplotype comparison is delegated to self.hap_comparator and must be at least
        self.threshold for a pair to be considered equal.
        
        The function emits tuples (i,j), where i and j are allele indices for which
        hapotypes id1-i and id2-j are equal.'''        
        return diff.equal_hap_pairs(self.__data, id1, id2, snp_start, snp_stop,
                                    self.__threshold, self.__hap_comparator)
        
####################################################################################
class IbdSegmentGlobalCacheIbdld:
    '''Maintains a cache of IBDLD results of IBD-1 sharing between genotyped samples.
    The results read from a file. The size of the cache is the number of pairs,
    which should be manageable as opposed to storing all segments.'''
    
    def __init__(self, file_name):
        '''Read an IBDLD lines from the file file_name and return a dictionary of
        (chrom,id1,id2)->1-based-line_number key-value pairs. Line numbers are to be
        subsequently used by linecache calls.'''
        with open(file_name, 'r') as f:
            reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
            line_number = {}
            count = 0
            for line in reader:
                if line:
                    count = count + 1
                    # Sort pair (i,j) so that i < j 
                    chrom = int(line[3])
                    id1 = int(min(line[1], line[2]))
                    id2 = int(max(line[1], line[2]))
                    line_number[(chrom, id1, id2)] = count
            self.file_name = file_name
            self.line_number = line_number

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def segments(self, chrom, id1, id2):
        '''Return the list of IBD segments shared by individuals id1, id2 on chrom chrom.'''
        a = min(id1, id2)
        b = max(id1, id2)
        try:
            line = linecache.getline(self.file_name, self.line_number[(chrom, a, b)]) 
            return np.array([segment for segment in self.__segments_in_line(line)])
        except KeyError:
            # No segments found
            return np.array([])
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __segments_in_line(self, line):
        '''A generator of segments from a given line.'''
        # Strip tab at end of line if exists
        items = line.rstrip('\n').split('\t')
        last_element = len(items) - 1 if items[-1] == '' else len(items)
        # Extract start [bp], stop [bp] for each segment
        for i in xrange(4, last_element, 6): yield [int(items[i]), int(items[i + 1])]

#---------------------------------------------
# Private Methods
#---------------------------------------------
def _hap_comparator_ibdld(a, b):
    '''Return a confidence measure in [0,1] of the equality of the haplotypes a and b.
    Defines the confidence of equality based on a model of IBDLD posterior probabilities.'''
    # TODO: replace this temporary implementation by ibdld_confidence()
    d = diff.hap_diff(a, b)
    return (1.0 * np.size(np.where(d == 0))) / max(1, np.size(np.where(d != INDETERMINATE)))
