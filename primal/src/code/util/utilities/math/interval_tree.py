'''
============================================================
Interval tree that allows fast queries of which intervals
intersect a location or an interval. ion among haplotypes).

Created on January 28, 2012
@author: http://hackmap.blogspot.com/2008/11/python-interval-tree.html
http://code.google.com/p/bpbio/source/browse/trunk/interval_tree/interval_tree.py
============================================================
'''
import operator

class IntervalTree(object):
    __slots__ = ('intervals', 'left', 'right', 'center', 'size')

    def __init__(self, intervals, depth=16, minbucket=64, _extent=None, maxbucket=512):
        """\
        `intervals` a list of intervals *with start and stop* attributes.
        `depth`     the depth of the tree
        `minbucket` if any node in the tree has fewer than minbucket
                    elements, make it a leaf node
        `maxbucket` even it at specified `depth`, if the number of intervals >
                    maxbucket, split the node, make the tree deeper.

        depth and minbucket usually do not need to be changed. if
        dealing with large numbers (> 1M) of intervals, the depth could
        be increased to 24.

        Usage:

         >>> ivals = [Interval(2, 3), Interval(1, 8), Interval(3, 6)]
         >>> tree = IntervalTree(ivals)
         >>> sorted(tree.find(1, 2))
         [Interval(1, 8), Interval(2, 3)]

        this provides an extreme and satisfying performance improvement
        over searching manually over all 3 elements in the list (like
        a sucker). 

        the IntervalTree class now also supports the iterator protocol
        so it's easy to loop over all elements in the tree:

         >>> import operator
         >>> sorted([iv for iv in tree], key=operator.attrgetter('start'))
         [Interval(1, 8), Interval(2, 3), Interval(3, 6)]


        NOTE: any object with start and stop attributes can be used
        in the incoming intervals list.
        """ 
        self.size = len(intervals)
        self.center = None
        depth -= 1
        if (depth == 0 or len(intervals) < minbucket) and len(intervals) < maxbucket:
            self.intervals = sorted(intervals, key=operator.attrgetter('start'))
            self.left = self.right = None
            return

        if _extent is None:
            # sorting the first time through allows it to get
            # better performance in searching later.
            intervals.sort(key=operator.attrgetter('start'))

        left, right = _extent or \
               (intervals[0].start, max(i.stop for i in intervals))
        # center = intervals[len(intervals)/ 2].stop
        center = (left + right) / 2.0

        self.intervals = []
        lefts, rights = [], []
        for interval in intervals:
            if interval.stop < center:
                lefts.append(interval)
            elif interval.start > center:
                rights.append(interval)
            else:  # overlapping.
                self.intervals.append(interval)

        self.left = lefts and IntervalTree(lefts, depth, minbucket, (intervals[0].start, center)) or None
        self.right = rights and IntervalTree(rights, depth, minbucket, (center, right)) or None
        self.center = center

    def find(self, start, stop):
        """find all elements between (or overlapping) start and stop"""
        if self.intervals and not stop < self.intervals[0].start:
            overlapping = [i for i in self.intervals if i.stop >= start and i.start <= stop]
        else:
            overlapping = []

        if self.left and start <= self.center:
            overlapping += self.left.find(start, stop)

        if self.right and stop >= self.center:
            overlapping += self.right.find(start, stop)

        return overlapping

    def __iter__(self):
        if self.left:
            for l in self.left: yield l

        for i in self.intervals: yield i

        if self.right:
            for r in self.right: yield r
   
    # methods to allow un/pickling (by pzs):
    def __getstate__(self):
        return { 'intervals' : self.intervals,
                'left'   : self.left,
                'right'  : self.right,
                'center' : self.center,
                'size'   : self.size }

    def __setstate__(self, state):
        for key, value in state.iteritems(): setattr(self, key, value)
    
class Interval(object):
    '''A closed interval [start,stop].'''
    __slots__ = ('start', 'stop')
    
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop
        
    def __repr__(self):
        return "Interval(%i, %i)" % (self.start, self.stop)

    def __key(self):
        '''Hash key.'''
        return (self.start, self.stop)

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
    
    def __getstate__(self):
        return {'start': self.start,
                'stop': self.stop }
    def __setstate__(self, state):
        for k, v in state.iteritems():
            setattr(self, k, v)

# if __name__ == '__main__':
#    import doctest
#    doctest.testmod()
