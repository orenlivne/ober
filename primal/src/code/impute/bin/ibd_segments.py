#!/usr/bin/env python
'''
============================================================
Run phasing in stages on a single chromosome part BED file. 

Created on August 16, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, impute as im, itertools, csv, optparse, traceback, util

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s [flags] <phased-data-file> <kinship-file>\n\n' \
        'Locate IBD segments among a subset of sample in an NPZ phased data set.\n' \
        'Sample pairs are read from standard input. Segments are written to standard output.\n' \
        '\tphased-data-file - NPZ file containing phasing results\n' \
        '\tkinship-file - Sorted identity coefficient file\n' \
        '\tpair-list-file - Sorted identity coefficient file\n' \
        '\tout-file - File to output segments to\n' \
        '\n' \
        'Example:\n' \
        'phased-data-file = chr22/hutt.phased.npz\n' \
        'kinship-file = hutt.kinship\n' \
        'pair-list-file contains the lines\n' \
        '0 1\n' \
        '...\n' \
        '0 100\n' \
        '\n' \
        'Type ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-v', '--debug', type='int', dest='debug', default=0,
                      help='Debug Level (0=quiet; 1=summary; 2=full debug)')
    parser.add_option('-i', '--input', type='str', dest='input_file', default=None,
                      help='Input file. If specified, reads from this file instead of from stdin')
    parser.add_option('-u', '--len-unit', type='str', dest='len_unit', default=im.PhaseParam().len_unit,
                      help='segment length to look for [cm|mbp]')
    parser.add_option('-l', '--min-len', type='str', dest='min_len', default=im.PhaseParam().min_len,
                      help='Minimum segment length to look for')
        
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 2:
        print usage
        sys.exit(1)
    phased_data_file, kinship_file = args

    try:
        # Load data
        problem = im.io.read_npz(phased_data_file)
        params = im.PhaseParam(kinship_file=kinship_file, debug=(options.debug >= 2), min_len=options.min_len)
        # Read all pairs from stdin first
        pairs = [(int(line[0]), int(line[1])) for line in 
                 csv.reader(open(options.input_file, 'rb') if options.input_file else sys.stdin,
                            delimiter=' ', skipinitialspace=True) if line]
        if options.debug >= 1:
            print params.min_len
            print params.debug
            print 'pairs', pairs
        
        # Loop over pairs and output segments to output file
        num_pairs = 4 * len(pairs)
        for k, ((i, j), (a, b)) in ((k, ((i, j), (a, b))) for (k, ((i, j), (a, b))) in 
                                    enumerate(itertools.product(pairs, itertools.product(im.constants.ALLELES, im.constants.ALLELES)))
                                    if (i, a) != (j, b)): # Note: "<" instead of "!=" should be implemented in pack_jobs.py predicates
            if options.debug >= 1:
                print 'Pair %d/%d: (%d,%d) (%d,%d)' % (k + 1, num_pairs, i, a, j, b)
            segments = im.ih.hap_segments_from_problem(problem, (i, a), (j, b), params)
            segments.save(sys.stdout)
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
