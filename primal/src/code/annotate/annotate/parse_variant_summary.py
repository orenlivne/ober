#!/usr/bin/env python
'''
============================================================
Parse variant summary mysql query result file.

Query: mysql -u hutt -phutt hutt -e "select is_qc, is_singleton, is_known, vartype, count(*) from hutt group by is_qc, is_singleton, is_known, vartype;"

Created on January 7, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import csv, sys

try:
    for x in csv.reader(sys.stdin, delimiter=','):
        print '%s_%s_%s_%s\t%s\t%s\t%s\t%s' % (x[21][3:], x[22], x[23], x[26], '1' if x[7] else '0', x[7] if x[7] else 0.0, x[24], x[25])
except (IOError, OSError):
    sys.exit(141)
