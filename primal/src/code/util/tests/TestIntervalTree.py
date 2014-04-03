'''
============================================================
Test interval trees.

Created on January 31, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, random, time, utilities.math.interval_tree as tr, cPickle, numpy as np
from numpy.testing.utils import assert_equal

class TestIntervalTree(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------       
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------       
    def test_query_time(self):
        '''Test interval query time with an interval tree vs. brute force.'''
        for scale in 2 ** np.arange(6, 11):
            self.__test_query_time_at_scale(scale)
    
    def test_serialization(self):
        '''Test that saving and loading a tree from a pickle file equals the original objectr.'''
        scale = 10000
        intervals, (START, STOP) = random_intervals(scale, 300)
        atree = tr.IntervalTree(intervals)
        btree = cPickle.loads(cPickle.dumps(atree, -1))
    
        af = sorted(atree.find(START, STOP)) 
        bf = sorted(btree.find(START, STOP))
        assert_equal(len(af), len(bf), 'Wrong # intervals in tree')
        for a, b in zip(af, bf):
            assert_equal(a.start, b.start)
            assert_equal(a.stop, b.stop)

    def test_query_is_correct(self):
        '''A use case in IBD segment sharing that seems to not have worked at some point.'''  
        intervals = [(36116080, 38818675), (16484792, 20993518)]
        START, STOP = (19056915, 29056915)
        tree = tr.IntervalTree(map(lambda x: tr.Interval(x[0], x[1]), intervals))
        res = tree.find(START, STOP)
        assert_equal(len(tree.intervals), 2, 'Wrong number of intervals in tree')
        print res
        assert_equal(res, [tr.Interval(16484792, 20993518)], 'Wrong query result')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __test_query_time_at_scale(self, scale):
        '''Test interval query time with an interval tree vs. brute force.'''
        intervals, (START, STOP) = random_intervals(scale) 
        
        t = time.time()
        bf = sorted(brute_force_find(intervals, START, STOP))
        btime = time.time() - t
    
        tree = tr.IntervalTree(intervals)
        tries = 10
        t = time.time()
        for i in xrange(tries):
            res = tree.find(START, STOP)
            assert_equal(sorted(res), bf)
        treetime = (time.time() - t) / tries
        self.assertTrue(treetime < 0.5 * btime, 'Query time not fast enough relative to brute-force')

        t = time.time()
        tree = tr.IntervalTree(intervals)
        for i in xrange(tries):
            bf = [i for i in intervals if i.stop >= START and i.start <= STOP]
        assert not set(bf).symmetric_difference(res) , (len(bf), len(res), set(bf).difference(res), START, STOP)
        assert_equal(sum(1 for _ in tree), len(intervals), 'iterator not working?')

#---------------------------------------------
# Methods
#---------------------------------------------
def random_intervals(scale, n=None):
    '''define a random set of intervals of size scale.'''
    if not n:
        n = 30 * scale
    intervals = [random_interval(scale) for _ in xrange(n)]
    intervals.append(tr.Interval(0, 50 * scale))
    return intervals, (39 * scale, 40 * scale)

def brute_force_find(intervals, start, stop):
    return [i for i in intervals if i.stop >= start and i.start <= stop]

def random_interval(scale):
    s = random.randint(1, 2000 * scale)
    return tr.Interval(s, s + random.randint(2, 6 * scale))
