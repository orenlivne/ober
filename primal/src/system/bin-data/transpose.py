#!/usr/bin/env python
'''Transpose matrix data.'''
import sys, csv, itertools
csv.writer(sys.stdout, delimiter=' ').writerows(itertools.izip(*csv.reader(sys.stdin, delimiter=' ')))
