#!/usr/bin/env python
'''
============================================================
Phase Hutterites in a single chromosome.

Created on August 16, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, impute as im, optparse, traceback, util
from impute.phasing.examples import OUT_PHASING

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
        'Phase Hutterites in a single chromosome.\n\n' \
        'Type ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-c', '--chr'          , type='int'           , dest='chrom', default=0,
                      help='Chromosome number to process. Overrides other options if non-zero. 0=do chromosomes in the range specified by the start, end chromosomes')
    parser.add_option('-s', '--start-chr'          , type='int', dest='start_chr', default=1,
                      help='Start Chromosome number')
    parser.add_option('-e', '--end-chr'          , type='int', dest='end_chr', default=22,
                      help='End Chromosome number')
    parser.add_option('-v', '--debug', type='int', dest='debug', default=0,
                      help='Debug Level (0=quiet; 1=summary; 2=full debug)')
    (options, args) = parser.parse_args(sys.argv[1:])    
    if len(args) != 0:
        print usage
        sys.exit(1)
    if (options.chrom < 0) or (options.chrom > 22):
        print usage
        print('\nMust specify a chromosome number in 1..22, or 0 for all chromosomes.')
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
        
    try:
        chromosomes = range(options.start_chr, options.end_chr + 1) \
        if options.chrom == 0 else [options.chrom]
        sys.stdout.write('Phasing chromosomes %s\n' % (repr(chromosomes,)))
        for chrom in chromosomes:
            sys.stdout.write('Phasing chromosome %d\n' % (chrom,))
            prefix = '%s/chr%d' % (OUT_PHASING, chrom)
            p = im.convert.main(pedigree=im.itu.HUTT_PED, prefix=prefix + '/hutt', npz=prefix + '/hutt.npz',
                                target='npz', debug=options.debug)
            im.phase.main(pedigree=im.itu.HUTT_PED, input=prefix + '/hutt.npz', output=prefix + '/hutt.phased.npz',
                          debug=options.debug, min_output=True)
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
