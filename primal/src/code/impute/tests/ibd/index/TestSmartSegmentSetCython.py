'''
============================================================
Test and time SmartSegmentSet implementation using native
Python vs. Cython.

Created on June 17, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import pyximport  # @UnresolvedImport
pyximport.install()

import time, unittest, test_util as tu, csv, itertools as it, impute as im
from impute.data.constants import MEGA_BASE_PAIR
from numpy.ma.testutils import assert_equal

class TestSmartSegmentSetCython(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    def setUp(self):
        '''Load segments from a test file.'''
        self.segments = map(lambda line: map(int, line),
                            csv.reader(open(tu.abs_path('segment/segments-0-100.out'), 'rb'),
                                       delimiter=' ', skipinitialspace=True))
        self.num_samples = 1415
        start, stop, min_len = 0, 99, 0.4 * MEGA_BASE_PAIR
        self.selected_segments = it.ifilter(lambda line: max(start, line[0]) < min(stop, line[1]) 
                                            and line[3] - line[2] >= min_len,
                                            self.segments)        
   
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_load_segments_native(self):
        '''Load segments from a segment list - native python implementation.'''
        start = time.time()
        s = im.smart_segment_set.SmartSegmentSet.from_list(self.num_samples, self.selected_segments)
        print time.time() - start, 'sec'
        assert_equal(s.size, 113660, err_msg='Wrong segment set size')

#     def test_load_segments_cython(self):
#         '''Load segments from a segment list - native python implementation.'''
#         start = time.time()
#         s = im.smart_segment_set_cython.SmartSegmentSet.from_list(self.num_samples, self.selected_segments)
#         print time.time() - start, 'sec'
#         assert_equal(s.size, 113660, err_msg='Wrong segment set size')

#---------------------------------------------
# Private Methods
#---------------------------------------------
