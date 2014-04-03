#!/usr/bin/env python
'''
============================================================
Convert Hutterite sample FINDIVs to our IDs. Reads from
stdin and converts each item.

Created on March 20, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, impute as im, os, numpy as np
from optparse import OptionParser

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <sample ID index file>\n' \
        'Convert CGI ID to sample ID.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--id-file', type='str'  , dest='id_file', default=os.environ['OBER'] + '/testdata/pedigree/hutterites.genotyped.tfam',
                      help='TFAM file of the genotyped Hutterites.')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 0:
        print usage
        sys.exit(1)
    
    # Read FINDIV-to-our-sample-index dictionary
    s = dict((v, k) for k, v in enumerate(np.loadtxt(options.id_file, usecols=[1], dtype=int)))
    
    # Convert lines
    try:
        for line in sys.stdin: print ' '.join(repr(s[int(x)]) for x in line.strip().split())
    except (IOError, OSError):
        sys.exit(141)
