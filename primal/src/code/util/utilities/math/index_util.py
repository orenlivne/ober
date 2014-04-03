'''
============================================================
Utility functions and classes.

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
# find_index
#---------------------------------------------
lib1 = ctypes.CDLL('%s/find_index.o.dylib' % DIR)
lib1.find_index.restype = ctypes.c_long
lib1.find_index.argtypes = (ctypes.c_long, ctypes.POINTER(ctypes.c_long), ctypes.c_long)

def first_index_long(a, value):
    '''Return the index of the first occurrence of a value in a dtype-long numpy array. If not found,
    returns -1.'''   
    return lib1.find_index(value, a.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), len(a))

#---------------------------------------------
# smallest_index
#---------------------------------------------
lib2 = ctypes.CDLL('%s/smallest_index.o.dylib' % DIR)

lib2.smallest_index_long.restype = ctypes.c_long
lib2.smallest_index_long.argtypes = (ctypes.c_long, ctypes.POINTER(ctypes.c_long), ctypes.c_long, ctypes.c_long, ctypes.c_long)

def first_occurrence_index_long(a, value, start, step):
    '''Return the index of the first occurrence of a value in a dtype-int numpy array. If not found,
    returns -1. If direction = 1, searches from the beginning of the array on; if direction=-1, searches
    from the last element back.'''
    return lib2.smallest_index_long(value, a.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), len(a), start, step)

lib2.smallest_index_byte.restype = ctypes.c_long
lib2.smallest_index_byte.argtypes = (ctypes.c_long, ctypes.POINTER(ctypes.c_int8), ctypes.c_long, ctypes.c_long, ctypes.c_long)

def first_occurrence_index_byte(a, value, start, step):
    '''Return the index of the first occurrence of a value in a dtype-int numpy array. If not found,
    returns -1. If direction = 1, searches from the beginning of the array on; if direction=-1, searches
    from the last element back.'''
    return lib2.smallest_index_byte(value, a.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)), len(a), start, step)
