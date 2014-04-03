#!/usr/bin/env python
'''
============================================================
Extract a set of lines from a file. The lines are read from
stdin. Output is written to stdout.

Created on October 4, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, sys

'''Loop over main file, select lines in the set ''lines'' and output them to the stream ''out''.'''
def write_file_subset(f, lines, out):
    for (k, line) in ((k, line) for k, line in enumerate(f) if k in lines):
        out.write(line)

'''Main program'''
if __name__ == '__main__':
    write_file_subset(open(sys.argv[1], 'rb'), set(np.loadtxt(sys.stdin, dtype=int)), sys.stdout)
