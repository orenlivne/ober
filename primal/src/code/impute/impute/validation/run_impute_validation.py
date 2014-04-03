#!/usr/bin/env python
'''
============================================================
Run an imputation validation experiment on a data set.
Data set can be anything that can be converted to a Problem
object, e.g., iplex (cf. imputation.reader).

Created on December 3, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, impute as im
from optparse import OptionParser
from impute.impute_test_util import HUTT_PED
from impute.phasing.examples import CHR22
 
####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <location_file> <true_type> <true_location> <save-location>\n' \
        'Impute a list of SNPs for which we have genotypes and compare with the imputed results.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-p', '--pedigree-file', type='str'  , dest='pedigree_file', default=HUTT_PED,
                      help='Pedigree file to read ALL subjects from (type+untyped, e.g., 3671 Hutterites)')
    parser.add_option('-g', '--genotype-id-file', type='str'  , dest='genotype_id_file',
                      default=CHR22 + '/hutt.tfam',
                      help='Pedigree file to read genotyped subject IDs from (e.g., 1415 Hutterites)')
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-s', '--snps', type='str', dest='snps', default=None,
                      help='Debug mode: comma-limited selected SNP indices to run')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 4:
        print usage
        sys.exit(1)
    # Parse custom SNP index list, if specified
    if options.snps: options.snps = [int(x) for x in options.snps.split(',')]
    location_file, true_type, true_location, save_location = args
    
    pedigree = im.io_pedigree.read(options.pedigree_file, options.genotype_id_file)
    p, t = im.v.iv.pipeline_validation_experiment(open(location_file, 'rb'), true_type, true_location, pedigree, debug=options.debug)
    stats = im.v.iv.stats_computer(t, p.g).stats()
    stats.summary()
    im.v.iv.plot_stats(stats, save_prefix=save_location + '/all' if save_location else None, fig_num=1, snp_style='discrete')
