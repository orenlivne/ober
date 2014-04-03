'''
============================================================
Low-level interval arithmetic and algorithms (related to
Identity-By-Descent (IBD) identification among haplotypes).

Created on August 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import util, numpy as np, networkx as nx, itertools, csv, itertools as it
from impute.data import constants
from impute.tools.genotype_tools import empty_errors_array
from types import GeneratorType
from collections import OrderedDict, Iterable
from blist import sortedlist
from impute.data.constants import ALLELES

#---------------------------------------------
# Constants
#---------------------------------------------
# Segment start and stop indices in a segment descriptor array [start,stop]
START, STOP, START_BP, STOP_BP = range(4) 

#---------------------------------------------
# Methods
#---------------------------------------------
def break_and_group_segments(segment_computer, debug=False):
    '''Return a generator of all (segment_snp_start,segment_snp_stop,h) contained in the snp
    region snp_window, where h is a subset of all haplotypes of the samples 'samples',
    whose elements share the segment.'''
    # Fetch IBD segments (each is a tuple (id1,id2,SNP range)
    segment_set = SegmentSet(segment_computer())
    segment_set.group_to_disjoint()
    return segment_set

def edges_to_segments(snps, edge, initial_phase, num_snps, output_phase=None, cover=False):
    '''Given a set of recombination locations in a sample (i.e., edges of some sort computed on a
    subset of the SNPs snps), return the set of recombination-free segments. 
    - edge is an index array into the snps array. 
    - An empty edge array return a single segment [0,num_snps]; 
    - if edge is None, returns an empty list.
    
    Options:
    If output_phase is not None, only output segments whose phase matches output_phase
    If cover = True, the IBD segments cover the entire 0..num_snps without holes; otherwise, they
    are restricted to cover the set snps.  
    '''
    # No recombination detected. Haplotype is entirely paternal or maternal
    #  ==> return a single segment
    if snps is None or edge is None:
        return []
    elif np.size(edge) == 0:
        # No recombination events
        segments = np.array([[0, num_snps, initial_phase]])
    else:
        # Recombination events exist. If deriv = -1 at an edge, it reflects a transition from non-equal
        # (1) to equal (0) paternal chromosome. So the parent's maternal chromosome is IBD with the
        # child; otherwise the paternal.
        # Fix first segment to start at the beginning of the chromosome rather than at the first het snp

        # Allocate IBD segment array
        num_edges = np.size(edge)
        segments = np.zeros((num_edges, 3), dtype=np.int)
        
        # First SNP starts at the beginning of the snps array. We cannot infer a segment outside the
        # input SNP set. 
        phase = initial_phase
        snp_left = 0 if cover else snps[0]
        for (count, edge_right) in enumerate(edge):
            segments[count, :] = [snp_left, snps[edge_right] + 1, phase]
            # Advance pointer to next het snp location, which might not be the nearest neighbor in the snps
            # array. There may be gaps of undetermined snps between two segments whose segment
            # membership we cannot ascertain.
            snp_left = snps[edge_right] + 1 if cover else snps[edge_right + 1]
            # Switch phase
            phase = 1 - phase

        # Last segment may be between the last edges and the last SNP in snps, which the segments are
        # are based on. We do not infer a segment outside the input SNP set if cover = False. 
        chrom_stop = num_snps - 1 if cover else snps[-1]  # and hence not num_snps-1, cf. above comment
        last_segment_exists = snps[edge[-1]] < chrom_stop
        if last_segment_exists:
            num_edges += 1
        segments = np.concatenate((segments, np.array([[snp_left, chrom_stop + 1, phase]])), axis=0)
    if output_phase:
        segments = segments[np.where(segments[:, 2] == output_phase)[0], :]
    return segments

def group_to_disjoint(segments, merge_consecutive=False):
    '''Given a -list- of key-value pairs ((start, stop), (v1,v2)), where the key denotes
    a segment [start,stop) of the infinite integer lattice on which v1 and v2 are equal,
    return a list of tuples ((a, b), U), where U is the union set of all values equal in [a,b).
    Here [a,b) are -disjoint- segments formed by intersections of the original segments,
    which might be partially-overlapping.
    
    If merge_consecutive is True, consecutive segments whose value sets (U for all U's)
    are identical are merged.'''

    sub_segments, intersections = form_intersections(segments)
    n = len(sub_segments)
    
    # For each bin, convert the pairwise adjacency list to an undirected graph.
    # Its connected components are the global hap sets that are equal in this segment.    
    if merge_consecutive:
        # Build a dictionary of (start, stop)-to-(set of connected components).
        disjoint_segments = OrderedDict(((sub_segments[i],
                                          frozenset(frozenset(component) for component 
                                                    in nx.connected_components(nx.Graph(intersections[i]))))
                                         for i in xrange(n)))
        
        # If two consecutive segments have identical sample sets, merge them
        if disjoint_segments:
            current_segment = None
            for (i, segment) in enumerate(disjoint_segments.iteritems()):
                if i == 0:
                    current_segment = segment
                else:
                    if segment[1] == current_segment[1]:
                        current_segment = ((current_segment[0][START], segment[0][STOP]), current_segment[1])
                        continue
                    else:
                        for component in current_segment[1]:
                            yield (current_segment[0], component)
                        current_segment = segment
            for component in current_segment[1]:
                yield (current_segment[0], component)
    else:
        # Yield the individual connected components of each segment without attempting to merge
        # consecutive segments
        for i in xrange(n):
            sub_segment = sub_segments[i]
            for component in nx.connected_components(nx.Graph(intersections[i])):
                yield (sub_segment, frozenset(component))


def primary_colors(segments, samples, num_primary):
    '''Given a -list- of key-value pairs ((start, stop), (v1,v2)), where the key denotes
    a segment [start,stop) of the infinite integer lattice on which v1 and v2 are equal,
    
    generate a haplotype coloring scheme from a list of IBD segments. samples is the set of 
    haplotype identifiers. if pair_gap is specified, haplotypes pairs are separated by this
    amount [pixels], assuming that their number is even. segment_gap is the # of SNPs to drop
    from each segment''s margins for the purpose of defining colors. This allows a slight
    overlap between segments without making them the same color, which is usually undesirable.
    
    Use a maximum of num_primary haps plus gray representing the rest.'''
    d = segments.to_group_to_color(samples=list(it.product(sorted(samples), ALLELES)))
    num_snps_in_color = np.array([sum(d[0][x][1] - d[0][x][0] for x, _ in y) for y in d[1]])
    ind = np.argsort(num_snps_in_color)
    return d[0], np.array(d[1][ind[-num_primary:]].tolist() + [[y for x in d[1][ind[:-num_primary]] for y in x]])

def group_to_color(segments, samples=None, segment_gap=25, snp_range=None):
    '''Given a -list- of key-value pairs ((start, stop), (v1,v2)), where the key denotes
    a segment [start,stop) of the infinite integer lattice on which v1 and v2 are equal,
    
    generate a haplotype coloring scheme from a list of IBD segments. samples is the set of 
    haplotype identifiers. if pair_gap is specified, haplotypes pairs are separated by this
    amount [pixels], assuming that their number is even. segment_gap is the # of SNPs to drop
    from each segment''s margins for the purpose of defining colors. This allows a slight
    overlap between segments without making them the same color, which is usually undesirable.'''
    
    # Define colors using slightly-smaller portions of the original segments
    segment_groups = __group_to_color([(__dilute_segment(x[0], segment_gap), x[1]) for x in segments], samples, snp_range=snp_range)

    # Translate back result to the original segments. Use sub-segments, dangling nodes of the
    # original segments (note: there may be more irrelevant dangling nodes in the diluted segments)
    sub_segments, _, _, segment_to_local_values, dangling_local_values = \
    form_intersections(segments, True, samples=samples, snp_range=snp_range)
    
    # Translate each connected component from a list of segments to a list of local values
    # (i.e., haplotype parts of the same color, for each color)
    return sub_segments, np.array([list(util.union_all(*(segment_to_local_values[x] for x in group))) 
                                   for group in segment_groups] + [[x] for x in dangling_local_values])

haps_of_color = lambda d, c: [(d[0][k], y[1]) for k, x in enumerate(d[1]) for y in x if y[0] == c]

def segments_with_value(A, value, L=0):
    '''Return a list of all maximal ranges [a,b) with b-a >= L such that A[a:b] = value. Useful for IBS
    segment identification. Note: A[0],A[-1] must not be nan.'''
    # Seemlessly handle boundary cases by appending values at both of A's ends,
    # forcing a change in A-value at the two boundaries, then subtracting one from the 'changes' indices
    # when calculating starts, ends
    changes = np.where(np.diff(np.concatenate(([np.nan], A, [np.nan])) == value))[0]
    starts, ends = changes[::2], changes[1::2]
    if L == 0:
        # Optimization
        return zip(starts, ends)
    else:
        lengths_minus_one = ends - starts
        idx = np.where(lengths_minus_one >= L)[0]  # ranges that are long enough
        return zip(starts[idx], ends[idx])

def flip_phase(segments):
    '''Flip the phase of all entries in a segment array.'''
    s = segments.copy()
    s[:, 2] = 1 - s[:, 2]
    return s

'''Return a logical array indicating whether the array values are in [x[0],x[1]) or not.'''
is_in_segment = lambda array, x: ((array >= x[0]) & (array < x[1]))

'''Return the sub-array of array whose values are in [x[0],x[1]).'''
in_segment = lambda array, x: array[np.where(((array >= x[0]) & (array < x[1])))[0]]

'''Flatten a list of tuples to a list of elements.'''
flatten = lambda x: [z for y in x for z in y]

'''Convert a SegmentSet s to a segment dictionary.'''
to_dict = lambda s: s.to_group_to_disjoint().to_group_by_snp_range()

'''Convert a list of segments obtained by SegmentSet.pprint_segments() back to a SegmentSet object.'''
segment_data_to_segment_set = lambda segment_data: SegmentSet([Segment(x[0], list(x[2]), (x[1][0], x[1][1])) for x in segment_data])

'''Return the stopping base pair location to be used in instantiating a Segment object. Attempts to
treat the chromosome right boundary case.'''
stop_bp = lambda bp, stop_snp, num_snps: bp[stop_snp] if stop_snp < num_snps else bp[num_snps - 1] + 1 

def segment_set_op(A, B, op):
    '''Union or intersect two sets of closed segments. op=''and' returns intersection, op='or' returns union.
    A set is specified by a list of tuples [(a[0],a[1]),...(a[2*n-2],a[2*n-1])] representing the segments
    [a[0],a[1]], ... [a[2*n-2],a[2*n-1]].'''
    # Operation to execute    
    if op == 'and': bit_op = lambda x, y: x & y
    elif op == 'or': bit_op = lambda x, y: x | y
    else: raise ValueError('Unsupported operation ''%s''' % (op,))
    
    a, b = flatten(A), flatten(B)

    # Create the list c, an encoded list of indices into the sorted union of a and b     
    i, j, n, m, c = 0, 0, len(a), len(b), []
    while i < n or j < m:
        if n > 0:
            bj = b[j] if j < m else a[-1] + 1  # Dummy value greater than all a-values if the b-list has been exhausted
            while i < n and a[i] <= bj:
                c.append(i)
                i += 1
        if m > 0:
            ai = a[i] if i < n else b[-1] + 1  # Dummy value greater than all b-values if the a-list has been exhausted
            while j < m and b[j] < ai:
                c.append(j + n)  # Encode b-values to be greater than all a-indices
                j += 1
    c, sz = np.array(c, dtype=int), len(c)
    
    # Identify the a-state: a_state[i] = 1 iff A contains the region c[i],c[i+1]. Similarly b.
    def to_state(ind):
        state = np.zeros((sz,), dtype=np.bool)
        for k in xrange(len(ind) / 2):
            s, t = ind[2 * k], ind[2 * k + 1]
            state[s:t] = 1
        return state
    a_state = to_state(np.where(c < n)[0])
    b_state = to_state(np.where(c >= n)[0])
        
    # Find ranges of consecutive 1's in the joint state array. Their endpoints are the desired output list
    decode = lambda x : a[x] if x < n else b[x - n]
    state = bit_op(a_state, b_state)
    return [(decode(c[x]), decode(c[y])) for x, y in segments_with_value(state, True)]

def union(A):
    '''Union of a set of closed segments.A set is specified by a list of tuples
    [(a[0],a[1]),...(a[2*n-2],a[2*n-1])] representing the segments [a[0],a[1]], ... [a[2*n-2],a[2*n-1]].'''
    a = list(A)
    if not a: return []
    (starts, ends), depth, union_endpoints = zip(*a), 0, []
    for point, point_type in sorted(it.chain(it.izip_longest(starts, [], fillvalue=START), it.izip_longest(ends, [], fillvalue=STOP))):
        if (depth == 0 and point_type == START) or (depth == 1 and point_type == STOP): union_endpoints.append(point)
        depth += (1 if point_type == START else -1)
    return zip(union_endpoints[0::2], union_endpoints[1::2]) 

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __dilute_segment(segment, segment_gap):
    '''Chop a margin of size segment_gap SNPs on either side of the segment tuple representing
    the SNP segment [a,b).'''
    start, stop = segment
    # Chop only if segment is not too small so that it would be empty after dilution
    return segment if stop - start <= 2 * segment_gap else (start + segment_gap, stop - segment_gap)

def __group_to_color(segments, samples=None, snp_range=None):
    '''Given a -list- of key-value pairs ((start, stop), (v1,v2)), where the key denotes
    a segment [start,stop) of the infinite integer lattice on which v1 and v2 are equal,
    Return a list of lists, each of contains segments of the same color (IBD sharing).'''
    (sub_segments, intersections, value_to_segments, _, _) = \
    form_intersections(segments, True, samples=samples, snp_range=snp_range)
    # Build a graph G where the nodes = segments and edges = (segments intersect AND their sample
    # sets intersect). G's connected components are groups, where group is a set of segments of the
    # same color.
    return nx.connected_components(nx.from_edgelist(itertools.chain.from_iterable(itertools.product(sharing_segments, sharing_segments)
              for sharing_segments in (util.union_all(*(value_to_segments[i][x] for x in component))
                                       for i in xrange(len(sub_segments))
                                       for component in nx.connected_components(nx.Graph(intersections[i]))))))
    
def form_intersections(segments, full_output=False, samples=None, snp_range=None):
    '''Return a list of bins whose ith element is the list of samples that are equal over
    [endpoints[i], endpoints[i+1]]. This is done by scanning all ibd_pairs entries
    and accumulating them into the intersecting segments' bins.
    Union set of all segment endpoints and sort ascending. These form
    the set of intersections of all pairwise_ibd segments.'''
    non_unique_endpoints = [s[0] for s in segments]
    # print 'snp_range', snp_range
    if snp_range is not None:
        non_unique_endpoints += [(snp_range[START], snp_range[STOP])]
    endpoints = sorted(reduce(set.union, non_unique_endpoints, set([])))
    
    # This index makes it easier to random-access bins
    endpoint_index = dict((x, index) for (index, x) in enumerate(endpoints))
    sub_segments = [(endpoints[i], endpoints[i + 1]) for i in xrange(len(endpoints) - 1)]
    # print 'endpoints', endpoints
    # print 'sub_segments', sub_segments
    
    n = len(endpoints)
    intersections = [list([]) for _ in xrange(n - 1)]
    if full_output:
        value_segment_tuples = [list([]) for _ in xrange(n - 1)]
        segment_local_value_tuples = []

    # Sweep segment intersections (bins) in order, accumulate segment info into each one
    for k, segment in enumerate(segments):
        snp = segment[0]
        for i in xrange(endpoint_index[snp[0]], endpoint_index[snp[1]]):
            intersections[i].append(segment[1])
            if full_output:
                value_segment_tuples[i] += [(x, k) for x in segment[1]]
                segment_local_value_tuples += [(k, (i, x)) for x in segment[1]]

    if full_output:
        # print 'segment_local_value_tuples', segment_local_value_tuples
        segment_to_local_values = util.to_set_dict(segment_local_value_tuples)
        # print 'segment_to_local_values', segment_to_local_values
        # print [x for s in segments for x in s[1]]
        all_samples = samples if samples is not None else util.union_all(x for s in segments for x in s[1])
        # print 'all_samples', all_samples
        all_local_values = set(itertools.product(xrange(n - 1), all_samples))
        dangling_local_values = all_local_values - util.union_all(*segment_to_local_values.values())
        # print 'dangling_local_values', dangling_local_values
        value_to_segments = [util.to_set_dict(d) for d in value_segment_tuples]
        # print 'value_to_segments', value_to_segments 
        util.to_set_dict((x, k) for k, s in enumerate(segments) for x in s[1])
        return (sub_segments, intersections,
                value_to_segments, segment_to_local_values, dangling_local_values)
    else:
        return (sub_segments, intersections)

####################################################################################
class Segment(object):  # implements interval_tree.Interval attributes
    '''Defines an IBD segment [start_snp,stop_snp) shared by set of samples or haplotypes. Note
    that this semi-closed-open interval.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, snp, samples, bp, error_snps=None, collapse_to_set=True, confidence=None,
                 cm=None):
        '''Initialize a segment [snp[0],snp[1]) that corresponds to the base pair segment
        [bp[0],bp[1]).

        
        If collapse_to_set=True, samples are collapsed into a set where each
        element appears only twice; otherwise, samples is stored intact.
        
        The optional array confidence of size snp[1]-snp[0] stores the probability that this
        segment has a certain property at each SNP (e.g., being IBD).

        cm = segment start and end in centi-Morgans.
        
        Note: the snp object can be mutable or immutable, dictating its subsequent mutability
        as a member of this object.'''
        if error_snps is None:
            error_snps = []
        self.__snp = snp
        self.__bp = bp
        if self.length == 0:
            raise ValueError('empty segment, [%d,%d), bp=(%d,%d)', snp[START], snp[STOP], bp[START], bp[STOP])
        self.samples = samples if not collapse_to_set or isinstance(samples, set) else set(samples)
        self.error_snps = error_snps
        self.confidence = confidence
        self.cm = cm

    #---------------------------------------------
    # Impl: interval_tree.Interval
    #---------------------------------------------
    @property
    def start(self):
        '''Return the segment (as Interval) starting point in base pairs.'''
        return self.__bp[START]
    
    @property
    def stop(self):
        '''Return the segment (as Interval) stopping point in base pairs.'''
        return self.__bp[STOP] - 1  # So that in terms of the Interval class' start, stop, the interval is closed: [start,stop]
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return 'Segment[[%4d,%4d], len=%6.2f, samples [%s]]' % \
            (self.__snp[START], self.__snp[STOP], self.length, ','.join(repr(s) for s in self.samples))
#        return 'Segment[[%4d,%4d], bp=[%d,%d], len=%6.2f, samples [%s]]' % \
#            (self.__snp[START], self.__snp[STOP], self.__bp[START], self.bp[STOP], self.length,
#             ','.join(repr(s) for s in self.samples))

    def __key(self):
        '''Hash key.'''
        return (self.__snp, self.samples)

    def __hash__(self):
        '''Hash key.'''
        return hash(self.__key())

    def __eq__(self, other):
        '''Equality of objects.'''
        return self.__key() == other.__key()

    def __ne__(self, other):
        '''Inequality of objects.'''
        return self.__key() != other.__key()

    def __cmp__(self, other):
        return cmp(self.__key(), other.__key())
    
    def middle_part(self, nearest_snp, bp, margin, collapse_to_set=True):
        '''Return a separate Segment instance where a margin is cut near both endpoints of the interval,
        leaving only the central portion of size (1-margin)*(original size) in base pairs. 
        If this portion is empty, returns None.
        
        round_position is a function of x that returns (snp, bp) where snp = the snp index nearest a
        bp position x, and bp is the corresponding base pair position of snp.'''
        m = 0.5 * margin
        bp_original = self.bp
        bp_middle = ((1 - m) * bp_original[0] + m * bp_original[1], m * bp_original[0] + (1 - m) * bp_original[1])
        (start_snp, stop_snp) = (nearest_snp(bp_middle[0]), nearest_snp(bp_middle[1]) + 1)
        if start_snp < stop_snp:
            # Non-empty segment left
            sub_segment = (start_snp, stop_snp)
            return Segment(sub_segment, self.samples, (bp[start_snp], bp[stop_snp]),
                           error_snps=in_segment(self.error_snps, sub_segment),
                           collapse_to_set=collapse_to_set)
        else:
            # Empty middle segment
            return None

    def sub_segment(self, sub_segment, bp, samples=None, collapse_to_set=True):
        '''Return a sub-segment [a,b).'''
        start_snp, stop_snp = sub_segment
        return Segment((start_snp, stop_snp), self.samples if samples is None else samples,
                        (bp[start_snp], bp[stop_snp]),
                        error_snps=in_segment(self.error_snps, sub_segment),
                        collapse_to_set=collapse_to_set)

    def remove_duplicate_samples(self):
        '''Remove samples that appear more than once in the sample set.'''
        self.samples = util.remove_duplicates(self.samples)

    def set_start(self, start_snp, bp_start):
        '''Set the starting SNP of the segment.'''
        # Tuples are immutable
        self.__snp = (start_snp, self.__snp[STOP])
        self.__bp = (bp_start, self.__bp[STOP])
        
    def set_stop(self, snp_stop, bp_stop):
        '''Set the stopping SNP of the segment.'''
        # Tuples are immutable
        self.__snp = (self.__snp[START], snp_stop)
        self.__bp = (self.__bp[START], bp_stop)
    
    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def snp(self):
        '''Return the SNP segment [start,stop] tuple.'''
        return self.__snp

    @property
    def bp(self):
        '''Return the segment endpoints in mega base pairs.'''
        return self.__bp

    @bp.setter
    def bp(self, bp):
        '''Set the bp position tuple.'''
        self.__bp = bp

    @property
    def num_snps(self):
        '''Return the segment length in SNPs.'''
        return self.__snp[STOP] - self.__snp[START]

    @property
    def length(self):
        '''Return the segment length in mega base pairs.'''
        return (self.__bp[STOP] - self.__bp[START]) / constants.MEGA_BASE_PAIR

    @property
    def length_cm(self):
        '''Return the segment length in centi-Morgans.'''
        return self.cm[STOP] - self.cm[START]

    @property
    def num_errors(self):
        '''Return the segment length in mega base pairs.'''
        return len(self.error_snps)

####################################################################################
class SegmentComposite(object):
    '''Defines an abstract set of IBD segments.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, segments=None):
        '''Initialize a segment set.'''
        # One underscore = protected member
        if segments is None:
            segments = []        
        self._segments = segments if isinstance(segments, list) else list(segments)
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------        
    def __iter__(self):
        '''Segment iterator.'''
        return self._segments.__iter__()
    
    def __getitem__(self, index):
        '''Get a segment in the list by index.'''
        return self._segments[index]

    def __len__(self):
        '''Get a segment in the list by index.'''
        return len(self._segments)

    def __radd__(self, other):
        '''Reverse-add a segment set to this object (binary operation).'''
        if isinstance(other, Iterable) and not other:
            # Empty + self = self
            return self
        elif isinstance(other, self.__class__):
            # SegmentComposite (of the same sub-class type as ours) is the only supported class
            return self +other
        else:
            raise ValueError('Cannot add an object of type SegmentComposite to a %s' % (other.__class__,))
        
    def __iadd__(self, segment):
        '''Append a segment or segment set to this object (unary operation). The self reference is not modified
        by this method.'''
        if isinstance(segment, self.__class__):
            # SegmentComposite (of the same sub-class type as ours)
            # self = self +segment # clean, but self is mutated
            self.append(segment)
        elif isinstance(segment, GeneratorType) or isinstance(segment, Iterable):
            # Collection of Segments
            for s in segment:
                self +=s
        else:#elif isinstance(segment, Segment):
            # Single segment
            self._segments.append(segment)
        #else:
            #raise ValueError('Cannot add an object of type %s to a SegmentComposite' % (segment.__class__,))
        return self

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property        
    def length(self):
        '''Return the number of items in this match set.'''
        return len(self._segments)

####################################################################################
class SegmentSet(SegmentComposite):
    '''Defines a set of IBD segments and associated genotype error locations.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, segments=None, errors=empty_errors_array()):
        '''Initialize a segment set.'''
        SegmentComposite.__init__(self, segments)
        self.errors = errors
        
    #---------------------------------------------
    # Operators
    #---------------------------------------------        
    def __repr__(self):
        return 'Segments:\n' + '\n'.join(repr(s) for s in self._segments) \
            + '\nErrors:\n' + '\n'.join(repr(e) for e in self.errors) if self._segments else ''
    
    def __add__(self, other):
        '''Merge two SegmentComposites (binary operation).'''
        if not other.errors.shape:
            errors = self.errors
        elif not self.errors.size:
            errors = other.errors
        else:
            errors = np.concatenate((self.errors, other.errors), axis=1)
        return SegmentSet(self._segments + other._segments, errors)
    
    def append(self, other):
        '''Same as __add__, only that the merged data is placed in self''s arrays.'''
        if not other.errors.shape:
            errors = self.errors
        elif not self.errors.size:
            errors = other.errors
        else:
            errors = np.concatenate((self.errors, other.errors), axis=1)
        self._segments = self._segments + other._segments
        self.errors = errors

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    @staticmethod        
    def load(inp):
        '''Load segments to the input stream inp. Works only when samples are haplotypes (i,a)
        where i=id and a=allele.'''
        return SegmentSet(Segment((items[0], items[1]),
                                  [(int(items[2 * i]), int(items[2 * i + 1])) for i in xrange(2, len(items) / 2)],
                                  (items[2], items[3]))
                          for items in ([int(x) for x in line] for line in csv.reader(inp, delimiter=' ')))  # @UndefinedVariable
    
    def save(self, out):
        '''Save segments to the output stream out.'''
        for s in self._segments:
            out.write('%d %d %d %d' % (s.snp[START], s.snp[STOP], s.bp[START], s.bp[STOP]))
            for sample in s.samples:
                out.write(' %s' % ' '.join(repr(x) for x in sample))
            out.write('\n')
            
    def pprint_segments(self, show_bp=False, show_cm=False):
        '''Return a pretty-printed-string of a list of segments.'''
        retval = '['
        n = self.length
        segments = self._segments
        for i in xrange(0, n):
            if (i > 0):
                retval += ' '
            s = segments[i]
            retval += '((%-4d, %-4d)' % (s.snp[START], s.snp[STOP])
            if show_bp:
                retval += ', (%-8d, %-8d, %7.3f, %d)' % (s.bp[START], s.bp[STOP], s.length, s.num_errors)
            if show_cm:
                retval += ', (%6.2f, %6.2f, %6.2f)' % (s.cm[START], s.cm[STOP], s.length_cm)
            retval += ', (%s))' % (','.join(repr(x) for x in s.samples) if len(s.samples) > 1 else '%s,' % repr(s.samples.__iter__().next())) 
            if (i < n - 1):
                retval += ',\n'
        retval += ']'
        return retval

    def to_sorted(self):
        '''Sort segments.'''
        return SegmentSet(sorted(Segment(x.snp, sorted(x.samples), x.bp, error_snps=x.error_snps,
                                         collapse_to_set=False, confidence=x.confidence) 
                                 for x in self._segments)) 

    def group_to_disjoint(self, merge_consecutive=False):
        '''Break segments (group by samples). Mutates the current instance.'''
        # Save BP positions of all SNPs in this set to be used in the broken segments' constructors
        snp_to_bp = dict(reduce(tuple.__add__, [((s.snp[START], s.bp[START]),
                                                 (s.snp[STOP], s.bp[STOP]))
                                                for s in self._segments], ())) 
        self._segments = list(Segment(snp, samples, (snp_to_bp[snp[START]], snp_to_bp[snp[STOP]]))
                              for (snp, samples) in 
                              group_to_disjoint([(s.snp, s.samples) for s in self._segments],
                                                merge_consecutive=merge_consecutive))
        return self

    def merge_consecutive(self):
        '''If two consecutive segments have identical sample sets, merge them.'''
        self._segments = list(self.__merge_consecutive())
        return self
        
    def __merge_consecutive(self):
        '''If two consecutive segments have identical sample sets, merge them.'''
        current_segment = None
        for (i, segment) in enumerate(self):
            if i == 0:
                current_segment = segment
            else:
                if segment.samples == current_segment.samples:
                    current_segment = Segment((current_segment.snp[START], segment.snp[STOP]),
                                              current_segment.samples,
                                              (current_segment.bp[START], segment.bp[STOP]),
                                              error_snps=current_segment.error_snps + segment.error_snps)
                    continue
                else:
                    yield current_segment
                    current_segment = segment
        yield current_segment

    def to_group_to_disjoint(self, merge_consecutive=False):
        '''Break segments (group by samples) in a new SegmentSet instance.'''
        # Save BP positions of all SNPs in this set to be used in the broken segments' constructors
        snp_to_bp = dict(reduce(tuple.__add__, [((s.snp[START], s.bp[START]), (s.snp[STOP], s.bp[STOP]))
                                                for s in self._segments], ()))
        segments = list(Segment(snp, samples, (snp_to_bp[snp[START]], snp_to_bp[snp[STOP]]))
                        for snp, samples in 
                        group_to_disjoint([(s.snp, s.samples) for s in self._segments],
                                          merge_consecutive=merge_consecutive))
        return SegmentSet(segments, errors=self.errors.copy())
    
    def to_group_to_color(self, samples=None, segment_gap=25, snp_range=None):
        '''Group to segment sets of unique color.'''
        return group_to_color([(s.snp, s.samples) for s in self._segments], samples=samples, segment_gap=segment_gap, snp_range=snp_range)
    
    def to_group_by_snp_range(self):
        '''Return a list of segment-SNP-range-to-list-of-sample-lists. Assumes that group_to_disjoint()
        has been called on this object. Segments and lists are lexicographically ordered.'''
        d = SnpRangeDictionary()
        for segment in self._segments:
            d.setdefault(segment.snp, sortedlist()).add(sortedlist(segment.samples))
        return d
        
    def remove_duplicate_samples(self):
        '''Remove samples that appear more than once in the sample set.'''
        for s in self.samples:
            s.remove_duplicate_samples()

    def segments_intersecting_snp(self, snp):
        '''Return the list of segments containing the SNP snp.'''
        return SegmentSet(x for x in self._segments if (x.snp[START] <= snp) and (snp < x.snp[STOP]))

    def segments_intersecting_bp(self, bp):
        '''Return the list of segments containing the base-pair position bp.'''
        return SegmentSet(x for x in self._segments if (x.bp[START] <= bp) and (bp < x.bp[STOP]))

    def segments_intersecting_bp_segment(self, bp):
        '''Return the list of segments containing the closed base-pair segment [bp[0], bp[1]].'''
        return SegmentSet(x for x in self._segments if max(x.bp[START], bp[START]) <= min(bp[STOP], x.bp[STOP] - 1))

    def segments_larger_than(self, min_length):
        '''Return the list of segments whose length >= min_length [base pairs].'''
        return SegmentSet(x for x in self._segments if x.length >= min_length)
    
####################################################################################
class SnpRangeDictionary(OrderedDict):
    '''Defines a dictionary: a map from a SNP range to a set of IBD sets.'''
    #---------------------------------------------
    # Operators
    #---------------------------------------------
    def __repr__(self):
        '''Return a pretty-printed-string of a list of segments.'''
        retval = '[\n'
        # n = self.length
        for (_, (snp, haps)) in enumerate(self.iteritems()):
            retval += '(%-4d, %-4d)\n' % (snp[START], snp[STOP])
            for ibd_set in haps:
                # retval += '\t\t[' + ','.join('(%-3d, %-1d)' % tuple(x) for x in ibd_set) + ']\n'
                retval += '\t\t' + repr(list(ibd_set)) + '\n'
        retval += ']'
        return retval

    #---------------------------------------------
    # Methods
    #---------------------------------------------

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def ibd2_segments(self):
        '''Return the list of segments for which one of the sharing samples is IBD-2 (homozygous).'''
        return [(k, x) for (k, v) in self.iteritems() for x in v if SnpRangeDictionary.__is_ibd2_segment(x)]

    @property
    def non_ibd2_segments(self):
        '''Segments that are not IBD-2 in any sharing sample.'''
        return [(k, x) for (k, v) in self.iteritems() for x in v if not SnpRangeDictionary.__is_ibd2_segment(x)]
    
    @property        
    def length(self):
        '''Return the number of items in this match set.'''
        return len(self)
    
    @property
    def all_samples(self):
        '''Return the union of all sample IDs from all sets of all segments.'''
        return set([]).union([hap for segment_list in self.itervalues() for s in segment_list for hap in s])

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    @staticmethod
    def __is_ibd2_segment(samples):
        '''Return true if and only if the segment set samples contains an IBD2 sample, i.e.,
        the haplotypes (x,0),(x,1) for some sample ID.'''
        return len(np.where(np.diff(sorted(y[0] for y in samples)) == 0)[0]) > 0
    
####################################################################################
class DisjointSegmentSet(object):
    '''Defines a light segment set that only holds the list of segment endpoints. Provides fast
    union, intersection operations. Segments are closed: [a,b], and must be non-overlapping.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, segments):
        '''Initialize from a list of segment endpoint 2-tuples.'''
        self._segments = segments
        
    #---------------------------------------------
    # Operators
    #---------------------------------------------
    def __repr__(self):
        '''Return a pretty-printed-string of a list of segments.'''
        return repr(self._segments)

    '''Hash and equality operators.'''
    def __key(self): return self._segments
    def __hash__(self): return hash(self.__key())
    def __eq__(self, other): return self.__key() == other.__key()
    def __ne__(self, other): return self.__key() != other.__key()

    '''Segment iterator.'''
    def __iter__(self): return self._segments.__iter__()
    '''Get a segment in the list by index.'''
    def __getitem__(self, index): return self._segments[index]

    def __radd__(self, other):
        '''Reverse-add a segment set to this object (binary operation).'''
        if not other:
            # Empty + self = self
            return self
        elif isinstance(other, self.__class__):
            # SegmentComposite (of the same sub-class type as ours) is the only supported class
            return self +other
        else:
            raise ValueError('Cannot add an object of type %s to a %s' % (self.__class__, other.__class__,))
            
    def __or__(self, other):
        '''Return the union of two segment sets (binary operation).'''
        return DisjointSegmentSet(segment_set_op(self._segments, other._segments, 'or'))

    def __and__(self, other):
        '''Return the intersection of two segment sets (binary operation).'''
        return DisjointSegmentSet(segment_set_op(self._segments, other._segments, 'and'))

    @property
    def size(self):
        '''Get the number of segments.'''
        return len(self._segments)
    
    @property
    def length(self):
        '''Get the total length of all segments. Can return long, not just int, unlike __len__.'''
        return sum((x[1] - x[0] + 1) for x in self._segments)
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
