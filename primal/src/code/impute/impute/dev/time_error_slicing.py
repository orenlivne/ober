'''
Created on Jan 17, 2013

@author: oren
'''
import numpy as np, timeit

def slice_1(a, rs, cs) :
    return a[rs][:, cs]

def slice_2(a, rs, cs) :
    return a[rs[:, None], cs]

rows, cols = 3218, 1415
a = np.random.rand(rows, cols) #@UnusedVariable
cs = np.array([2, 1]) #@UnusedVariable
rs = np.arange(rows) #@UnusedVariable

calls = 10
print timeit.timeit('slice_1(a, rs, cs)',
                    'from __main__ import slice_1, a, rs, cs',
                    number=calls) / calls
print timeit.timeit('slice_2(a, rs, cs)',
                    'from __main__ import slice_2, a, rs, cs',
                    number=calls) / calls
