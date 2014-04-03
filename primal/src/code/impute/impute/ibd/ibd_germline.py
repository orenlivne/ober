'''
============================================================
IBD segment identification between haplotypes using
the GERMLINE algorithm.

Created on September 14, 2012
@author: Oren Livne <livne@uchicago.edu>
@see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2652213/pdf/318.pdf
============================================================
'''
import itertools, numpy as np, util
from impute.tools import genotype_tools as gt
from impute.data import constants
from hashable import Hashable
from impute.ibd.segment import Segment, START, STOP, START_BP, SegmentSet
from types import GeneratorType
from impute.tools.param import PhaseParam
from impute.tools.pedigree_tools import parent_type_str

# TODO: replace with a HET snp threshold in PhaseParam
FILL_THRESHOLD = 0.95
# TODO: replace with a size dependent on the #meioses between the most distant relatives we're considering
DEFAULT_WORD_SIZE = 100

#---------------------------------------------
# Methods
#---------------------------------------------
####################################################################################
def ibd_germline(problem, samples):
    '''Return all IBD segments among the genotyped samples' haplotypes using GERMLINE. Segments are
    rounded to the nearest slice.'''
    ibd_computer = GermlineIbdComputer(PhaseParam())  
    h_mat = _HapMatrix(problem, gt.genotyped(problem, samples))
    m = ibd_computer.ibd_segments(h_mat)
    m.group_to_disjoint()
    return m

class GermlineIbdComputer:
    '''Computes IBD segments among a haplotype set using the GERMLINE algorithm.''' 

    def __init__(self, params=PhaseParam()):
        self.params = params
        self.debug = params.debug
        # TODO: move into params
        self.min_segment_length = 0.0
        # TODO: calculate min_segment_length from s_ibd, s, slice_size
        self.max_difference = 2
        
    def ibd_segments(self, h):
        '''GERMLINE Algorithm 4:
        Yield all IBD segments between samples in a haplotype matrix.
        
        initial_slice_size I.
        * Chance of obtaining a false segment in the initial slice among N samples: 
                            1-(1-2**-I)**(N*(N-1)/2) (e.g., I=20, N=10: 4.3e-5)
        * Chance of missing a segment in the initial slice among N samples with genotype error rate e: 
                            1-(1-e)**(2*I*N) (e.g., I=20, N=10, e=0.01: .86. But at least a -short- segment is missing only.)
        
        min_segment_length L - minimum segment length [base pair].
        '''        
        # Divide h into slices of size slice_size abutting at SNP indices e
        if self.debug:
            h.print_info()
        h.form_slices(self.params.initial_slice_size, self.params.slice_size)
        if self.debug:
            print 'Sliced'
            print h.endpoints
        
        # Find segments in first slice
        k = 0
        (m, mprev_extended) = (SegmentSet([]), h.to_match_set(k, self._match(h.get_slice(k))))
        if self.debug:
            print 'Initial segments'
            print mprev_extended

        # Loop over remaining slices 
        m_extended = mprev_extended
        for k in xrange(1, h.num_slices + 1):
            if k < h.num_slices:
                if self.debug:
                    print '=' * 40, 'Slice', k, '=' * 40
                word = h.get_slice(k)
                mk = self._match(word)
                if self.debug:
                    print 'Slice', k, 'Range', h.get_slice_endpoints(k), \
                    'SNPs', h.get_slice_endpoint_snps(k), 'Matches', len(mk), mk
                    print 'Previous matches'
                    print mprev_extended
                m_extended = self.__extend(h, k, mprev_extended, mk)
                if self.debug:
                    print 'Extended matches'
                    print m_extended
                self.__extend_partial_forward(word, mprev_extended, m_extended, self.max_difference)
                if self.debug:
                    print 'Extended-partial matches'
                    print m_extended
            # Add terminated matches to final match set
            if self.debug:
                print 'Adding terminal matches to final match set:'
                print mprev_extended
            # Segment objects are assumed to be semi-closed-open: [a,b)
            m += (Segment((h.original_snps[h.snps[s.snp[START]]], h.original_snps[h.snps[s.snp[STOP]]] + 1),
                          set([h.index_to_hap(s.samples[0]), h.index_to_hap(s.samples[1])]),
                          s.bp, collapse_to_set=False)
                  for s in mprev_extended.itervalues() if s.length >= self.min_segment_length)
            mprev_extended = m_extended
        if self.debug:
            print 'Final match set:'
            print m.pprint_segments(True)
        return m
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------    
    def _match(self, h):
        '''GERMLINE Algorithm 1:
        Return a set of all matches (i,j) of row indices for which the matrix h is identical.'''
        # Create a dictionary of haplotypes-to-(their-row-indices) "bins" 
        d = {}
        for (i, hi) in enumerate(h):
            d.setdefault(Hashable(hi), []).append(i)
        if self.debug:
            print 'Initial matches dict:'
            for (k, v) in d.iteritems():
                print k.unwrap().astype(np.int), v
            
        # Union all pairs (i,j) from each bin into the returned set.
        return reduce(set.union,
                      (set(((i, j) for (i, j) in itertools.product(row_set, row_set) if i < j))
                      for row_set in d.itervalues()))
    
    def __extend(self, h, k, mprev_extended, m):
        '''GERMLINE Algorithm 2:
        Extend matches m_extended from previous slice with the matches m of the current slice.
        Extended matches are moved from m_extended to the returned segment set, so this function
        mutates m_extended.'''
        m_extended = h.to_match_set(k, m)
        if self.debug:
            print 'New matches'
            print m_extended
        for (pair, match) in [(p, m) for (p, m) in m_extended.items() if mprev_extended.has_key(p)]:
            existing = mprev_extended[pair]
            if self.debug:
                print 'Extending match', match, 'to start at', existing.snp[START], existing.bp[START]
            match.set_start(existing.snp[START], existing.bp[START])
            mprev_extended -= existing
        return m_extended
    
    def __extend_partial_forward(self, word, m_prev, m_curr, max_difference):
        '''GERMLINE Algorithm 3:
        Augment the remaining matches in the previous slice m_prev with partial matches with up to
        max_difference mismatches in the m_curr slice.'''
        (stop, stop_bp) = (m_curr.slice_snp[STOP], m_curr.slice_bp[STOP])
        # print m_prev.items()
        for match in [m for (p, m) in m_prev.items() 
                      if util.num_differences(word[p[0], :], word[p[1], :]) <= max_difference]:
            if self.debug:
                print 'Extending partial match', match, 'to stop at', stop, stop_bp
                p = m.samples
                print '# diffs', p, util.num_differences(word[p[0], :], word[p[1], :])
            match.set_stop(stop, stop_bp)
            m_prev -= match
            m_curr += match

####################################################################################
class _HapMatrix(object):
    '''Extract the haplotypes of the n samples in the array 'sample' and stack them as a
    binary 2n x s matrix H, where rows 2*i and 2*i+1 correspond to individual sample[i].
    True means allele = 1, False means allele = 2 in the original haplotype.
    
    If parent_type is specified, H is n x s and consists only of parent_type-type haplotypes.
    
    NOTE: this object maintains a copy (not a reference of the original haplotype data
    array h. It clears missing data from all rows and returns the remaining column indices as
    the second arguemnt in the returned tuple.'''
    '''Holds haplotype IBD segment information.'''    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, problem, sample, snps=None, parent_type=None):
        '''Construct a haplotype matrix for a sample array sample.'''
        self.__problem = problem
        self.__sample = sample
        self.__parent_type = parent_type
        h = problem.haplotype.data
        # Original SNPs = s
        s = snps if snps is not None else problem.snp_range
        self.original_snps = s
        H = (h[:, sample, :].reshape((h.shape[0], 2 * len(sample))).transpose() if parent_type is None \
        else h[:, sample, parent_type].transpose())[:, s]
        # Restricted SNP set where all samples have hap data
        self.snps = np.where(np.min(H, axis=0) != constants.MISSING)[0]
        # Corresponding BP positions
        self.__bp = problem.haplotype.snp['base_pair'][s[self.snps]]
        # Haplotype binary matrix
        self.h = (H[:, self.snps] == 1)
        # Convert self.h row index s to original haplotype index (sample_index, hap_parent_type) 
        self.index_to_hap = lambda t: (self.__sample[t / 2], t % 2) if parent_type is None \
        else (self.__sample[t], parent_type)

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return 'HapMatrix[%d x %d, sample=%s, sample_index=%s]' % \
            (self.shape[0], self.shape[1], repr(self.__sample),
             repr([self.__problem.pedigree.sample_index[x] for x in self.__sample]))
    
    def print_info(self):
        '''Print descriptive information.'''
        print 'Hap Matrix (%d x %d), samples %s, parent_type %s' % \
            (self.shape[0], self.shape[1], repr(tuple(self.sample)), parent_type_str(self.__parent_type))
        print 'Hap matrix row-to-sample mapping:'
        for i in xrange(self.shape[0]):
            print '%d -> (%d,%d)' % ((i,) + self.index_to_hap(i))

    def form_slices(self, initial_slice_size, slice_size):
        '''Mark slices in H.'''
        initial_slice_size = min(self.num_snps, initial_slice_size)
        slice_size = min(self.num_snps, slice_size)
        boundary = [0] + range(initial_slice_size, self.num_snps, slice_size)
        # If last slice is too short, merge it with its predecessor
        if self.num_snps - boundary[-1] <= 0.7 * slice_size:
            boundary[-1] = self.num_snps
        else:
            boundary += [self.num_snps]
        endpoints = np.array([(boundary[i], boundary[i + 1] - 1) for i in xrange(0, len(boundary) - 1)])
        self.__endpoints = np.concatenate((endpoints, self.__bp[endpoints]), axis=1)

    def get_slice(self, k):
        (start, stop) = self.__endpoints[k][START:START_BP]
        return self.h[:, start:stop]

    def get_slice_endpoints(self, k):
        return self.__endpoints[k]

    def get_slice_endpoint_snps(self, k):
        '''Convert to original SNP indices.'''
        return (self.snps[self.__endpoints[k][0]], self.snps[self.__endpoints[k][1]])

    def to_match_set(self, k, matches):
        return _MatchSet(self.__endpoints[k][START:START_BP], self.__endpoints[k][START_BP:],
                         (self.__new_segment(k, pair) for pair in matches))
        
    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def endpoints(self):
        return self.__endpoints

    @property
    def sample(self):
        '''Return the sample index corresponding to hap rows.'''
        return self.__sample
        
    @property
    def shape(self):
        '''Return the entire hap matrix's shape.'''
        return self.h.shape

    @property
    def num_snps(self):
        '''Return the number of SNPs (# columns).'''
        return len(self.snps)

    @property
    def num_slices(self):
        return self.__endpoints.shape[0]
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __new_segment(self, k, pair):
        '''Return a Segment instance for slice k and the sample pair 'pair'.'''
        # Create mutable copies of SNP endpoints, bp
        return Segment(tuple(self.__endpoints[k][START:START_BP]), pair,
                       tuple(self.__endpoints[k][START_BP:]), collapse_to_set=False)

####################################################################################
class _MatchSet(object):
    '''Defines a set of IBD segment matches.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, slice_snp, slice_bp, segments):
        '''Initialize a segment set covering up to slice slice_stop.'''
        self.slice_snp = tuple(slice_snp)
        self.slice_bp = tuple(slice_bp)
        self.__dict = dict(((s.samples, s) for s in segments)) 
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------        
    def __repr__(self):
        '''Pretty-print a list of segments.'''
        n = self.length
        if n == 0:
            return '[]'
        result = '%d matches:\n[' % (n,)
        for (i, s) in enumerate(self.__dict.itervalues()):
            if (i > 0):
                result += ' '
            result += '((%-4d, %-4d), %s, (%-8d, %-8d), %7.3f)' % \
            (s.snp[START], s.snp[STOP], repr(s.samples), s.bp[START], s.bp[STOP], s.length)
            if (i < n - 1):
                result += ',\n'
        result += ']\n'
        return result

    def has_key(self, pair):
        '''Check if pair is in the segment list.'''
        return self.__dict.has_key(pair)

    def __getitem__(self, pair):
        return self.__dict[pair]
    
    def itervalues(self):
        return self.__dict.itervalues()

    def items(self):
        return self.__dict.items()
    
    def __add__(self, segment):
        '''Add segment.'''
        if isinstance(segment, GeneratorType):
            for s in segment:
                self = self +s
        else:
            self.__dict[segment.samples] = segment
        return self

    def __iadd__(self, segment):
        return self +segment
    
    def __isub__(self, segment):
        del(self.__dict[segment.samples])
        return self
    
    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property        
    def length(self):
        '''Return the number of items in this match set.'''
        return len(self.__dict)
    
#---------------------------------------------
# Private Methods
#---------------------------------------------
def _filled_samples(problem, sample, threshold=FILL_THRESHOLD):
    '''Return the family members with filled haplotypes.''' 
    fill = problem.fill_fraction(sample=sample)
    return np.intersect1d(fill[fill[:, 1] > threshold, 0].astype('int'), sample)

def _filled_members(problem, family, threshold=FILL_THRESHOLD):
    '''Return the family members with filled haplotypes.''' 
    return _filled_samples(problem, gt.genotyped_members(problem, family), threshold)
