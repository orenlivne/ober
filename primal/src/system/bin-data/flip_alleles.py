#!/usr/bin/env python
'''
============================================================
Flip 1 and 2 alleles in a PLINK input TPED for a subset
of the SNPs.

Created on May 24, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv, db_gene
from optparse import OptionParser

# Flip allele number. Missing values are kept intact.
recode = lambda x: 0 if x == 0 else 3-x

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s\n' \
        'Flip 1 and 2 allele coding in a TPED stdin.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-l', '--snps', type='str'  , dest='snps', default=None,
                      help='If specified, this is the file name containing a subset of selected SNPs to flip. If not specified, all SNPs are flipped.')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 0:
        print usage
        sys.exit(1)

    # Read SNP list
    snps = set(map(lambda x: x.rstrip('\r\n').rstrip('\n'), open(options.snps, 'rb').readlines())) if options.snps else None

    for line in csv.reader(sys.stdin, skipinitialspace=True, delimiter=' '):
        if not snps or line[1] in snps:
            print ' '.join(line[:4]) + ' ' + ' '.join(map(str, map(recode, map(int, line[4:]))))
        else:
            print ' '.join(line)
