#!/usr/bin/env python
'''
============================================================
Extract CEU columns for the annotation database from the
annovar annotation fields of the 96 Hutterites variants.

Reads from stdin, writes to stdout.

Output format:

hash_code is_ceu maf_ceu ref_allele alt_allele xref var_region var_func var_mutation
22_16050035_16050036_snp        0       0.0     A       C       -       intergenic      NONE(dist=NONE),POTE\
H(dist=206296) MISSENSE
22_16050407_16050408_snp        1       0.06    T       C       rs2844883       intergenic      NONE(dist=NO\
NE),POTEH(dist=205924) NONSENSE
...

Created on December 18, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import csv, sys

try:
    for x in csv.reader(sys.stdin, delimiter=','):
        var_mutation = x[34] if len(x) > 34 and x[34] else '-'
        var_category = x[0] if x[0] else '-'
        print '%s_%s_%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % \
            (x[21][3:],
             x[22],
             x[23],
             x[26],
             '1' if x[7] else '0',
             x[7] if x[7] else 0.0,
             x[24],
             x[25],
             x[8] if x[8] else '-',
             var_category,
             x[1] if x[1] else '-',
             var_mutation,
             x[14] if x[14] else 0,
             x[15] if x[15] else '-')
except (IOError, OSError):
    sys.exit(141)
