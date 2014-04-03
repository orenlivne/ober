#!/usr/bin/env python
'''
============================================================
Impute SNPs from 98 CGI whole-genome-sequencing samples to
all Hutt samples using IBD cliques calculated and indexed
from Affymetrix SNP data.

Created on February 4, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, impute as im
from optparse import OptionParser

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <phased-problem.npz> <ibd-segment-file> <out-file.npz>\n' \
        'Impute all (or some) SNPs in the phased data set. Save results in NPZ format.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-s', '--snps', type='str', dest='snps', default=None,
                      help='Debug mode: comma-limited selected SNP indices to run')
    parser.add_option('-i', '--index-file', type='str'  , dest='index_file',
                      default=os.environ['OBER'] + '/data/cgi/README.assembly_sample_subject.csv',
                      help='CGI+FINDIV ID index file')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 3:
        print usage
        sys.exit(1)
    # Parse selected SNP index list
    if options.snps: options.snps = [int(x) for x in options.snps.split(',')]

    # Run imputation        
    _, t = im.v.iv.impute_problem(args[0], args[1], options.index_file, snps=options.snps, debug=options.debug)
    # Save results
    t.save(args[2])
    
