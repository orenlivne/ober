'''
============================================================
IBD segment identification between a proband and distant
relatives in the pedigree graph.

Created on August 12, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, impute as im
from impute.data.constants import PATERNAL, MATERNAL, ALLELES
from impute.ibd import ibd
from types import GeneratorType
from impute.ibd.segment import SegmentSet, START, STOP, Segment
from impute.tools.pedigree_tools import surrogate_parents
from impute.tools.param import PhaseParam
from impute.ibd.ibd import genotype_ibs_segments
from utilities.math.index_util import first_occurrence_index_byte

#---------------------------------------------
# Methods
#---------------------------------------------
def prob_ibd_ibs(problem, id1, id2, snps, params):
    '''Is the segment s (defined by the SNP array snps, which may be a frame within the actual segment)
    an IBD segment between samples id1 and id2 or not? Outputs an IBD probability estimate.

    This estimate is based on IBS, i.e., the probability is ~1.0 at all SNPs within IBS segments
    (decaying a bit near the boundaries), and 0.0 outside.'''
    ibs_segments = genotype_ibs_segments(problem.genotype, id1, id2, snps,
                                         min_ibs_len_snp=params.min_ibs_len_snp,
                                         debug=params.debug)
    cm = problem.info.snp['dist_cm']
    p_ibd = np.zeros_like(snps, dtype=np.bool)
    for segment in ibs_segments:
        start, stop = segment.snp
        p_ibd[start:stop] = _ibs_confidence(cm[start:stop], segment.length)
    return p_ibd

def ibd_segments_with_surrogate_parents(problem, i, min_path_length, max_path_length,
                                        surrogate_parent_fill_threshold=0.9,
                                        params=PhaseParam(),
                                        prob_ibd_calculator=prob_ibd_ibs):
    '''A utility method that calculates likely IBD segments between i and surrogate parents.
    Supports minimum required IBD segment length vs. d meioses separating two samples. Assumes
    # an exponential distribution and computes the value x for which P(x >= X) = ibd_length_upper_percentile.
     
    Important: use a short median filter window size (~3) to catch more spikes since we are
    less certain that spikes are genotype errors than in parent-child IBD calculations
    '''
    # Find the set J of phased relatives of m=max_path_length proximity; if not found, increment
    # m to at most max_path_length until such are found
    # relatives = self.filled_relatives(i, self.max_path_length, self.het_fill_threshold)
    relatives = RelativeCollection.in_neighborhood(i, problem, min_path_length, max_path_length,
                                                   params.surrogate_parent_fill_threshold)
    return ibd_segments_with_relatives(problem, i, relatives.info['index'], params, prob_ibd_calculator)

def ibd_segments_with_relatives(problem, i, relatives, params, prob_ibd_calculator, use_kinship=True):
    '''Calculates likely IBD segments between i and a list of relatives.'''
    computer = IbdDistantSegmentComputer(problem, i, params, prob_ibd_calculator)
    return computer.segments(relatives, use_kinship=use_kinship)

def best_segment_index(p, segments):
    '''Return an array with the segment index within a SegmentSet segments
    to use for phasing at each SNP. Requires the SNP base-pair location array.'''
    # bp = p.info.snp['base_pair']
    # Holds the highest segment ranking index at each SNP
    rank = np.zeros(p.num_snps, dtype=np.float) 
    # Holds segment index to use for phasing at each SNP
    s = -np.ones(p.num_snps, dtype=np.int)
    # Loop over segments; update ranks within to this segment if they are higher than the highest
    # rank so far
    for i, segment in enumerate(segments):
        start, stop = segment.snp
        # L = segment.length
        # x = bp[start:stop]
        # phi =  L * ibd_ld_confidence((x - 0.5 * (x[0] + x[-1])) / MEGA_BASE_PAIR, L)
        phi = segment.confidence
        to_update = np.where(phi >= rank[start:stop])
        # print 'Updating segment to %d at %s' % (i, repr(to_update))
        rank[start:stop] = np.maximum(phi, rank[start:stop])
        s[start:stop][to_update] = i
    return s

####################################################################################
class IbdDistantSegmentComputer(object):
    '''Computes segments between sample i (phased or unphased) and multiple phased surrogate parents j.'''
    
    def __init__(self, problem, i, params, prob_ibd_calculator=prob_ibd_ibs):
        '''Initialize a segment computer for sample i.'''
        self.problem = problem
        self.i = i
        self.params = params
        self.prob_ibd_calculator = prob_ibd_calculator
        self.num_snps = problem.num_snps
        
        # Setup that depends only on i
        # Find sufficiently-long IBD segments between i and j, for each j in J. The length threshold
        # is decreased like the expected segment length in base pairs, i.e., O(1/m) where
        # m=shortest path length between i and j.
        g, h = problem.data
        self.het = h[im.gt.where_heterozygous(g, i), i, :]
        # Fetch frames of the chromosome in question. Assuming all SNPs are on the same chromosome.
        self.frames = problem.frames[problem.info.snp['chrom'][0]]
        self.bp = problem.info.snp['base_pair']        
        self.cm = problem.info.snp['dist_cm']
        self.gi_exists = im.recode.filled_genotype(g[:, i])

    def segments(self, relatives, use_kinship=True):
        '''Yield IBD segments between a sample i and a relative j. relatives is an Nx2 array of
        [relative_index, relatedness_path_len_to_i].
         
        Supports minimum required IBD segment length vs. d meioses separating two samples. Assumes
        an exponential distribution and computes the value x for which P(x >= X) = ibd_length_upper_percentile.
        
        Important: use a short median filter window size (3-5) to catch more spikes since we are
        less certain that spikes are genotype errors than in parent-child IBD calculations
        
        If is_i_phased = True, will try to match i's haplotypes against j's. Otherwise, will match
        one of i's haps against j's genotype (assuming both i's haps are equal, and phased only at hom SNPs).  
        '''
        # Check i's het fill only once per call, not once per entry in the relatives list 
        is_i_phased = (self.het_fill_fraction > self.params.het_fill_threshold)
        if np.isscalar(relatives): return self.__segments(relatives, is_i_phased)
        else: return sum([self.__segments(j, is_i_phased, use_kinship=use_kinship) for j in relatives], SegmentSet([]))
        
    #---------------------------------------------
    # Operators
    #---------------------------------------------
    @property
    def het_fill_fraction(self):
        '''Current fill fraction of i at het SNPs.'''     
        return 1.0 * np.size(np.where(self.het != im.constants.MISSING)[0]) / max(1, np.size(self.het))

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __segments(self, j, is_i_phased, use_kinship=True):
        '''Yield IBD segments between a sample i and a phased, POO-determined relative j.
         
        Supports minimum required IBD segment length vs. d meioses separating two samples. Assumes
        an exponential distribution and computes the value x for which P(x >= X) = ibd_length_upper_percentile.
        
        Important: use a short median filter window size (3-5) to catch more spikes since we are
        less certain that spikes are genotype errors than in parent-child IBD calculations
        
        If is_i_phased = True, will try to match i's haplotypes against j's. Otherwise, will match
        one of i's haps against j's genotype (assuming both i's haps are equal, and phased only at hom SNPs).  
        '''
        g, h = self.problem.data
        is_i_phased = (self.het_fill_fraction > self.params.het_fill_threshold)
        
        # Kinships between i's parents and j's parents for parent-of-origin delineation in the produced
        # segments (we need to know which of i's haplotypes is IBD with the j haplotype that fits an IBD=1
        # segment)
        sample_id = self.problem.pedigree.sample_id
        if use_kinship:
            predecessors = self.problem.pedigree.graph.predecessors
            i_parents_ids = [sample_id[predecessors(self.i)[a]] for a in ALLELES]
            j_parents_ids = [sample_id[predecessors(j)[a]] for a in ALLELES]
            K = np.array([[self.params.kinship(x, y) for y in i_parents_ids] for x in j_parents_ids])
        if self.params.debug:            
            print '-' * 70
            print 'IBD Segments between i and a surrogate parent j'
            print 'i=%d ID %8d fill %.3f%% Phased? %s' % \
            (self.i, sample_id[self.i], 100.*self.problem.fill_fraction_of_sample(self.i), "yes" if is_i_phased else "no") 
            print 'j=%d ID %8d fill %.3f%%' % (j, sample_id[j], 100.*self.problem.fill_fraction_of_sample(j))
            u, d = im.pt.lowest_common_ancestor(self.problem.pedigree.graph, self.i, j)
            if u:
                print 'Lowest-common ancestor: %d, depth=%d' % (u, d)
            print '-' * 70
            print 'Frame lengths', [len(x) for x in self.frames]

        #--------------------------------------------
        # Pick a frame of independent SNPs 
        #--------------------------------------------
        # SNPs with full Gi, Gj data
        mask_gi_gj_exist = self.gi_exists & im.recode.filled_genotype(g[:, j])
        gi_gj_exist = np.where(mask_gi_gj_exist)[0]
        # Find Largest-intersecting independent SNP frame
        pair_frames = [np.intersect1d(frame, gi_gj_exist) for frame in self.frames]
        frame_number = np.argmax(np.array([len(x) for x in pair_frames]))
        frame_size = len(frame)
        # Add the first and last SNP in the domain that have i,j data for full frame coverage
        start = first_occurrence_index_byte(mask_gi_gj_exist, True, 0, 1)
        stop = first_occurrence_index_byte(mask_gi_gj_exist, True, len(gi_gj_exist) - 1, -1) 
        frame = np.sort(np.union1d([start, stop], pair_frames[frame_number]))
        if self.params.debug:
            print 'Frame #', frame_number, ', size', frame_size, 'start', start, 'stop', stop, frame.tolist()
        
        #--------------------------------------------
        # Calculate IBD posterior and IBD mask 
        #--------------------------------------------
        prob_ibd = self.prob_ibd_calculator(self.problem, self.i, j, frame, self.params)
        is_ibd = prob_ibd > 0.0
        if self.params.debug:
            # print 'IBD mask', is_ibd.astype(np.uint)
            print 'Raw IBD segments', [(frame[x[START]], frame[x[STOP]] if x[STOP] < len(frame) else self.num_snps)
                                       for x in ibd.segment.segments_with_value(is_ibd, True, 3)] 
            
        #--------------------------------------------
        # Calculate Haplotype matching mask 
        #--------------------------------------------
        # Interpolate prob_ibd from frame SNPs to all SNPs (0 outside IBD segments; linear interpolation
        # within segments). Interpolate mask_is_ibd (piecewise constant: 1 within segments, 0 outside).
        mask_prob_ibd = np.zeros((self.num_snps,), dtype=np.float)
        mask_is_ibd = np.zeros((self.num_snps,), dtype=np.bool)
        for fr, full_segment in ((fr, self.__frame_segment_to_snp(frame, fr))
                                 for fr in ibd.segment.segments_with_value(is_ibd, True)):
            start, stop = full_segment[START], full_segment[STOP]  # Original indices of segment boundaries
            fr_start, fr_stop = fr[START], fr[STOP]  # Frame indices of segment boundaries
            mask_prob_ibd[start:stop] = np.interp(self.cm[start:stop], self.cm[frame[fr_start:fr_stop]], prob_ibd[fr_start:fr_stop])
            mask_is_ibd[start:stop] = True
        
        # i is phased ==> calculate difference mask between each i-hap and j-hap
        # i is not phased ==> calculate difference between i's genotype and each j-hap. Output a dummy i-hap in segments
        mask_is_hap_fit = ((im.diff.all_diffs(h, self.i, j) if is_i_phased else \
        np.array([im.recode.ibs_diff_gh(g[:, self.i], h[:, j, j_allele]) for j_allele in ALLELES])) <= 0)
        mask = mask_is_ibd & mask_is_hap_fit
        
        # filtered_mask = ibd.filter_diff(mask, 'median', 5)
        # Filter segments based on length
        ibd_segments = SegmentSet([])
        # If use_kinship=True and i is unphased, find which of i's parents has higher kinship to the
        # j_allele-parent of j. Assign the segment to be between the corresponding i and j haplotypes.
        #
        # If use_kinship=False, assign the IBD segment to an arbitary haplotype (PATERNAL). Inferior
        # and should only be used during debugging. 
        for k, (i_allele, j_allele) in enumerate([(i_allele, j_allele) 
                                                  for i_allele in (ALLELES if is_i_phased else [np.argmax(K[:, j_allele]) if use_kinship else PATERNAL]) 
                                                  for j_allele in ALLELES]):
            mask_pair = mask[k]
            filtered_mask = ibd.filter_diff(mask_pair, 'median', 5)
            error = (filtered_mask != mask_pair)
            # Output only segments of 1's of at least a trivially-small length
            # print filtered_mask
            segments = ibd.segment.segments_with_value(filtered_mask, True, 3)  
            if self.params.debug:
                print '--- IBD segments between (%d,%d), (%d,%d)' % (self.i, i_allele, j, j_allele)
            min_len = self.params.min_len  # @UnusedVariable
            len_unit = self.params.len_unit
            long_segments = [Segment(s, [(self.i, i_allele), (j, j_allele)], (bp_start, bp_stop),
                                     error_snps=np.where(ibd.segment.is_in_segment(error, s)),
                                     confidence=mask_prob_ibd[s[START]:s[STOP]],
                                     cm=(cm_start, cm_stop), collapse_to_set=False)
                             for s, bp_start, bp_stop, cm_start, cm_stop in
                             ((t,
                               self.bp[t[START]], im.segment.stop_bp(self.bp, t[STOP], self.num_snps),
                               self.cm[t[START]], im.segment.stop_bp(self.cm, t[STOP], self.num_snps)
                               ) for t in segments)
                             if ((cm_stop - cm_start >= min_len) if (len_unit == 'cm') else (bp_stop - bp_start >= im.constants.MEGA_BASE_PAIR * min_len))]
            ibd_segments += long_segments
            if self.params.debug:
                print 'Long segments (>= %.2f %s):' % (min_len, len_unit)
                print '\n'.join('\t' + repr(segment) for segment in long_segments) if long_segments else '\t-'
                print ''
        return ibd_segments

    def __frame_segment_to_snp(self, frame, segment):
        '''Convert a segment=[a,b) in frame coordinates to SNP indices of the original problem.'''
        return (frame[segment[START]], frame[segment[STOP]] if segment[STOP] < len(frame) else self.num_snps)

#---------------------------------------------
# Private Methods
#---------------------------------------------
def _ibs_confidence(x, L):
    '''A rough estimated confidence of IBD within an IBS segment whose endpoints are
    in centi-Morgans. x = array of base-pair positions of segment SNPs.
    
    Note: confidence values range from 0 to infinity, not from 0 to 1.'''
    return L * im.ibdld.ibd_ld.ibd_ld_confidence((x - 0.5 * (x[0] + x[-1])), L)

####################################################################################
class RelativeCollection:
    '''A helper class that holds information about a set of relatives of a sample.'''
    
    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, sample, info=None):
        '''Initialize a relative collection of the proband ''sample''.'''
        self.sample = sample
        self.info = info if info else \
        np.empty((0,), dtype=[
                        ('index', np.uint),  # Sample index = unique identifier within a Genotype object
                        ('distance', '(2,)i4')  # Minimum # meioses between relative and proband through father, mother
                        ])

    @staticmethod
    def from_data(sample, relatives_data):
        '''Create a RelativeCollection from a list of full sibs of sample.'''
        r = RelativeCollection(sample)
        r += [x for x in relatives_data]
        return r

    @staticmethod
    def from_sibs(sample, sibs):
        '''Create a RelativeCollection from a list of full sibs of sample.'''
        r = RelativeCollection(sample)
        r += [(x, (2, 2)) for x in sibs]
        return r
    
    @staticmethod
    def in_neighborhood(sample, problem, min_path_length, max_path_length, fill_threshold):
        '''Create a RelativeCollection with genotype relatives filled to fill_threshold of
        at most d=max_path_length proximity.'''
        genotyped_relatives = dict((x, y) for (x, y) in 
                                   surrogate_parents(problem.pedigree.graph, sample,
                                                     min_depth=min_path_length,
                                                     max_depth=max_path_length).iteritems()
                                   if problem.is_genotyped(x))
        f = problem.fill_fraction(sample=genotyped_relatives.keys())
        r = RelativeCollection(sample)
        if f.size:
            filled = dict(zip(f[:, 0].astype('int'), f[:, 1]))
            r += [(x, y) for (x, y) in genotyped_relatives.iteritems()
                  if filled[x] >= fill_threshold]
            # Sort by ascending priority (increasing depth, then decreasing fill %)
            # return a[np.lexsort((a[:, 1], -a[:, 2]))] if a.size else None
        return r
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def min_segment_length(self, ibd_length_upper_percentile):
        '''Return a list of (index, l) tuples, where l = minimum expected IBD segment length vs.
        d meioses separating the proband. Assumes an exponential segment length distribution L,
        and computes the value l for which P(l >= L) = ibd_length_upper_percentile.'''
        return zip(self.info['index'], -100.0 * np.log(ibd_length_upper_percentile) / self.min_distance) 
    
    #---------------------------------------------
    # Operators
    #---------------------------------------------
    def __add__(self, relative):
        '''Add relative info (either a single tuple (index, (distance_paternal, distance_maternal)),
        or a list/generator expression of such tuples).'''        
        if isinstance(relative, GeneratorType):
            relative = list(relative)
        current_size = self.info.size
        if isinstance(relative, list): 
            # Bulk insert
            self.info.resize(current_size + len(relative))
            self.info[current_size:] = relative
        else:
            # Single-record insert
            self.info.resize(current_size + 1)
            self.info[-1] = relative
        return self

    def __iadd__(self, segment):
        return self +segment
    
    def __repr__(self):
        return 'RelativeCollection[' + repr(self.info) + ']'
     
    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def length(self):
        '''# of relatives.'''
        return len(self.info)
    
    @property
    def min_distance(self):
        '''Return the minimum of the paternal and maternal pedigree distances between each relative
        and the proband.'''
        return np.minimum(self.info['distance'][:, PATERNAL], self.info['distance'][:, MATERNAL])

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
