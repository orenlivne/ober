'''
============================================================
Digital filters.

To compile (Unix):
make

Created on Jun 1, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import ctypes, os

#---------------------------------------------
# Constants
#---------------------------------------------

# Current directory, for Cython library absolute paths
DIR = os.path.dirname(__file__)

#---------------------------------------------
# filter
#---------------------------------------------

lib1 = ctypes.CDLL('%s/filter.o.dylib' % DIR)
lib1.median_filter.restype = ctypes.c_long
lib1.median_filter.argtypes = (ctypes.POINTER(ctypes.c_long), ctypes.c_long, ctypes.c_long, ctypes.c_long)

def median_filter(a, filter_size):
    '''Apply a median filter.'''   
    return lib1.median_filter(a.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), a.shape[0], a.shape[1], filter_size)
