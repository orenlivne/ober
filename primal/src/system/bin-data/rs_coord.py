#!/usr/bin/env python
'''
============================================================
Return the hg19 chromosome and coordinate of an RS number.
If no arguments are specified, reads from stdin, otherwise
reads from file. One SNP per line.

Created on February 20, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import re, sys, itertools as it
from Bio import Entrez

Entrez.email = 'god@heaven.gom'
ENTREZ_CHUNK_SIZE = 1000  # Max number of records to pass to a batch Entrez query 
COORD_REGEX = re.compile('assembly=GRCh37.p10 \| chr=(\d+) \| chr-pos=(\d+)')

def grouper(iterable, n, fillvalue=None):
    '''Collect data into fixed-length chunks or blocks.
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx'''
    args = [iter(iterable)] * n
    return it.izip_longest(fillvalue=fillvalue, *args)

####################################################################################
if __name__ == '__main__':
    input_file = open(sys.argv[1], 'rb') if len(sys.argv) == 2 else sys.stdin
    for chunk in grouper(input_file, ENTREZ_CHUNK_SIZE):
        term = ','.join(x.strip() for x in it.takewhile(lambda x: x, chunk))
        response = Entrez.efetch(db='SNP', id=term, rettype='flt', retmode='flt').read().split('\n')
        for chr, pos in (m.groups() for line in response for m in COORD_REGEX.finditer(line)): print chr, pos
