'''
============================================================
Test IBD segment queries using interval trees.

Created on January 31, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, impute as im, tempfile
from numpy.testing.utils import assert_equal

class TestSmartSegmentSet(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------       
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------       
    def test_interval_query(self):
        '''Test interval queries.'''
        s = im.smart_segment_set.SmartSegmentSet.from_list(4, TestSmartSegmentSet._segment_data())
        assert_equal(len(s), 4, 'Wrong # entries in set')
        assert_equal(s.size, 8, 'Wrong # segments in set')
        assert_equal(TestSmartSegmentSet.segment_ids(s.find(20, 50, (0, 0))), set([1, 2, 3, 5, 6]), 'Wrong intersecting segment set')
        assert_equal(TestSmartSegmentSet.segment_ids(s.find(20, 50, (0, 0), [(1, 0), (1, 1), (2, 0)])), set([1, 2, 3, 5]), 'Wrong intersecting segment set')
        assert_equal(TestSmartSegmentSet.segment_ids(s.find(20, 50, (0, 0), [(1, 0), (1, 1), (2, 1)])), set([1, 2, 3]), 'Wrong intersecting segment set')
        assert_equal(TestSmartSegmentSet.segment_ids(s.find(20, 50, (1, 1), [(1, 1), (2, 1)])), set([7]), 'Wrong intersecting segment set')
        
    def test_serialization(self):
        '''Test that saving and loading a tree from a pickle file equals the original objectr.'''
        s = im.smart_segment_set.SmartSegmentSet.from_list(4, TestSmartSegmentSet._segment_data())
        out_file = tempfile.TemporaryFile()
        s.save(out_file)
        out_file.seek(0)
        s2 = im.smart_segment_set.SmartSegmentSet.load(4, out_file)
        out_file.close()
        assert_equal(s, s2, 'Saving and loading did not restore the original the sample set')
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    @staticmethod
    def _segment_data():
        '''Return a list of test segments. Row format: snp_range bp_range hap1 hap2'''
        # snp_start is used here as a unique segment identifier
        return [[0, 1, 0, 10, 0, 0, 1, 1],
                [1, 1, 40, 100, 0, 0, 1, 1],
                [2, 1, 5, 30, 0, 0, 1, 0],
                [3, 1, 40, 60, 0, 0, 1, 0],
                [4, 1, 70, 100, 0, 0, 1, 0],
                [5, 1, 0, 100, 2, 0, 0, 0],
                [6, 1, 10, 60, 3, 1, 0, 0],
                [7, 1, 10, 70, 1, 1, 2, 1]]

    @staticmethod
    def segment_ids(segments):
        '''snp_start is used here as a unique segment identifier.'''
        return set(map(lambda x: x.snp_start, segments))
