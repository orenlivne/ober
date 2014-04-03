#!/usr/bin/env python
'''
============================================================
Run phasing in stages on a single chromosome part BED file. 

Created on August 16, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, util
from impute import impute_test_util as itu
from optparse import OptionParser
from impute.preprocess import convert
from impute.phasing import phase

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <in-dir> <out-dir> <chr> <part>\n' \
        'Run phasing on an individual chromosome part BED file.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-s', '--stage', default=0, type='int',
                      help='Run only this phasing stage (0=all; 1=trivial; 2=family: parent-offspring; \
                      3=family: children comparison; 4: family: sib comparison);')
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    (options, args) = parser.parse_args(sys.argv[1:])
    if len(args) != 4:
        print usage
        sys.exit(1)

    # in_dir =  'phasing.20121130/split_chr'
    # out_dir = 'phasing.20121130/individual'
    # chrom = 5
    # part = 0
    in_dir, out_dir, chrom, part = args[0], args[1], int(args[2]), int(args[3])
    print 'Running phasing in stages'
    print 'in_dir  = %s' % (in_dir,)
    print 'out_dir = %s' % (out_dir,)
    print 'chrom   = %d' % (chrom,)
    print 'part    = %d' % (part,)

    out_dir = '%s/chr%d' % (out_dir, chrom)
    util.mkdir_if_not_exists(out_dir)
    
    npz_file = '%s/hutt.stage0.npz' % (out_dir,)
    if not os.path.exists(npz_file) and options.stage == 0:
        convert.main(pedigree=itu.HUTT_PED,
                     prefix='%s/chr%d/hutt_chr%d_part%d' % (in_dir, chrom, chrom, part),
                     npz=npz_file, target='npz', debug=True)
    
    for stage in (range(1, 5) if options.stage == 0 else [options.stage]):
        phase.main(pedigree=itu.HUTT_PED,
                          input='%s/hutt.stage%d.npz' % (out_dir, stage - 1),
                          output='%s/hutt.stage%d.npz' % (out_dir, stage),
                          stage=stage, debug=options.debug)
