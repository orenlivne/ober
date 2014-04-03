#!/usr/bin/env python
'''
============================================================
Run phasing in stages on a single chromosome part BED file. 

Created on August 16, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, impute as im, csv, optparse, traceback, util

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s [flags]\n\n' \
        'Test python import.\n' \
        '\n' \
        'Type ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-i', '--input', type='str', dest='input_file', default=None,
                      help='Input file. If specified, reads from this file instead of from stdin')
        
    (options, args) = parser.parse_args(sys.argv[1:])
    if len(args) != 0:
        print usage
        sys.exit(1)

    try:
        print 'Start'
        im.examples.test()
        for line in (line for line in csv.reader(open(options.input_file, 'rb') if options.input_file else sys.stdin,
                                                 delimiter=' ', skipinitialspace=True) if line):
            print line
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
