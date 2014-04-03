# encoding: utf-8
# filename: smart_segment_set_cython.py
'''
============================================================
Defines a set of IBD segments between sample pairs.
Optimized for access speed. Native python implementation.

Created on August 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
from utilities.math import interval_tree
from impute.data.constants import MEGA_BASE_PAIR

#---------------------------------------------
# Methods
#---------------------------------------------
'''Return an interval tree set for a key. Standardize each segment x to have the hap ''key''
listed first in the sample list.'''
_interval_set = lambda key, segments: \
    interval_tree.IntervalTree([x if x.sample0 == key else  
                                PairSegment(x.sample1, x.sample0, x.snp_start, x.snp_stop, x.bp_start, x.bp_stop)
                                for x in segments])

####################################################################################
class PairSegment(interval_tree.Interval):
    __slots__ = ('sample0', 'sample1', 'snp_start', 'snp_stop', 'bp_start', 'bp_stop')
    
    '''Defines an IBD segment [start_snp,stop_snp) shared by a pair of haplotypes. Note
    that this semi-closed-open interval.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, sample0, sample1, snp_start, snp_stop, bp_start, bp_stop):
        '''Initialize a segment [snp[0],snp[1]) that corresponds to the base pair segment
        [bp[0],bp[1]).'''
        # Note: the Interval class represents a closed interval [bp_start,bp_stop]
        interval_tree.Interval.__init__(self, bp_start, bp_stop - 1)
        self.sample0, self.sample1 = sample0, sample1
        self.snp_start, self.snp_stop = snp_start, snp_stop
        self.bp_start, self.bp_stop = bp_start, bp_stop
        if bp_stop <= bp_start:
            raise ValueError('empty segment, snp=[%d,%d), bp=(%d,%d)', snp_start, snp_stop, bp_start, bp_stop)

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self): return 'PairSegment[[%4d,%4d], len=%6.2f, samples [%s,%s]]' % \
                (self.snp_start, self.snp_stop, len(self), self.sample0, self.sample1)
    '''Object equality.'''
    def __key(self): return (self.__snp, self.samples)
    def __hash__(self): return hash(self.__key())
    def __eq__(self, other): return self.__key() == other.__key()
    def __ne__(self, other): return self.__key() != other.__key()
    def __cmp__(self, other): return cmp(self.__key(), other.__key())
    '''Return the segment length in mega-base pairs.'''
    def __len__(self): return (self.bp_stop - self.bp_start) / MEGA_BASE_PAIR

    '''Methods to allow un/pickling.'''
    def __getstate__(self):
        return {'sample0': self.sample0, 'sample1': self.sample1,
                'snp_start': self.snp_start, 'snp_stop': self.snp_stop,
                'bp_start': self.bp_start, 'bp_stop': self.bp_stop }

    def __setstate__(self, state):
        for key, value in state.iteritems(): setattr(self, key, value)
    
####################################################################################
class SmartSegmentSet(object):
    '''Defines a segment set suitable for fast intersection queries. Implemented using haplotype
    indexing. Each segment is inserted twice, under the keys of both of its sharing haplotypes.
    Each key''s segment list can be efficiently searched using an interval tree.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, n, segments):
        '''Initialize a segment set for n samples numbered 0..n-1 from a PairSegment list.'''
        # PairSegment index: a 1-D array with entries corresponding to (0,0),(0,1),...,(n-1,0),(n-1,1),
        # where n=#samples. Each entry is a list of segments between the hap (i,j) and all other haps.

        # Accumulate segments into lists
        self._n = n
        self._index = [[] for _ in xrange(2 * n)]
        for x in segments:
            self(x.sample0).append(x)
            self(x.sample1).append(x)
        # Convert segment lists to searchable objects
        self._index = [_interval_set((k / 2, k % 2), x) for k, x in enumerate(self._index)]
        
    #---------------------------------------------
    # Operators
    #---------------------------------------------
    '''Return a pretty-printed-string of a list of segments.'''
    def __repr__(self): return '[\n' + '\n'.join('(%4d,%4d): %s' % (k / 2, k % 2, repr(x)) for k, x in enumerate(self._index)) + '\n]'
    '''Object equality.'''
    def __key(self): return self._index
    def __hash__(self): return hash(self.__key())
    def __eq__(self, other): return self.__key() == other.__key()
    def __ne__(self, other): return self.__key() != other.__key()
    '''Get the number of segments in the set.'''
    def __len__(self): return self._n
    '''Get an index entry corresponding to a haplotype key=(sample,allele).'''    
    def __call__(self, key): return self._index[np.uint(2 * key[0] + key[1])]
    '''Get the number of entries = # haplotypes in the set. Recall that every segment appears twice.'''
    @property
    def size(self): return sum(x.size for x in self._index) / 2

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    @staticmethod       
    def from_list(n, s):
        '''Load segments from a segment data list.
        Each segment is assumed to be between a PAIR of haplotypes only.
        Only list-type ''samples'' fields are supported.'''
        return SmartSegmentSet(n, [PairSegment((x[4], x[5]), (x[6], x[7]), x[0], x[1], x[2], x[3]) for x in s]) 

    @staticmethod       
    def load(n, inp):
        '''Load segments from the file inp. 
        Each segment is assumed to be between a PAIR of haplotypes only.
        Only list-type ''samples'' fields are supported.'''
        return SmartSegmentSet.from_list(n, np.loadtxt(inp, dtype=np.uint))

    def save(self, out, samples=None):
        '''Save segments to the output stream out. 
        Each segment is assumed to be between a PAIR of haplotypes only.
        Only list-type ''samples'' fields are supported.'''
        s = np.array([(x.snp_start, x.snp_stop, x.bp_start, x.bp_stop) + x.sample0 + x.sample1 
                      for lst in self._index for x in lst._segments
                      if x.samples[0] < x.samples[1]])
        np.savetxt(out, s, fmt='%d')
        
    def among_samples(self, samples):
        '''Return a smaller segment set instance with segments among a subset of the samples.'''
        return SmartSegmentSet(max(samples) + 1,
                               [x for lst in self._index for x in lst._segments
                                if x.samples[0] < x.samples[1]
                                and x.samples[0][0] in samples and x.samples[1][0] in samples])

    def of_haps(self, haps):
        '''Return a smaller segment set instance with segments between haps and all other haps.'''
        return SmartSegmentSet(self._n, filter(lambda x: x.samples[0] < x.samples[1],
                                               sum((self(x)._segments for x in haps), [])))

    def find(self, start, stop, source, target=None):
        '''find all segments or overlapping the segment [start,stop] (units=base-pairs)
        between the haplotype ''source'' and haplotypes in ''target''. If target=None, searches in all
        haplotypes.'''
        intersecting = self(source).find(start, stop)
        return [s for s in intersecting if s.samples[0] == source and s.samples[1] in target] if target is not None \
                else intersecting

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
