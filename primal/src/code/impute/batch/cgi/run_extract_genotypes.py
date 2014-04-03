#!/usr/bin/env python
'''
============================================================
Convert Hutterites WGS data at SNP locations read from
standard input to PLINK TPED format.

Prerequisites: tabix

Created on November 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, util, traceback, impute as im
from optparse import OptionParser
from impute.cgi.extract_genotypes import extract_genotypes

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s [options] <input_file> <output_file>\n' \
        'Extract SNP data from tabixed CGI var files and convert to PLINK recoded-12 TPED or\n' \
        'python NPZ format. Base-pair ranges are read from input_file in the format: chr start_bp end_bp\n' \
        '* Multiple ranges can be specified (one per line).\n' \
        '* If input_file = ''-'', input is taken from standard input.\n' \
        '* Use 23 for chrX, 24 for chrY. Chrs XY, M are not yet supported.\n' \
        '* By default, allele labels are swapped so that 1 is always the major and 2 is always the minor.\n' \
        '\n' \
        'Examples:\n' \
        'Extract part of chr X to a PLINK TPED:\n\techo "23 0 1000001" | %s -t plink - hutt-x-cgi\n' \
        'Extract part chr 1 to an npz file, don\'t swap major-minor allele labels:\n\techo "1 0 1000001" | cgi2plink.py -n -t plink - hutt-1-cgi\n' \
        '\n' \
        'Type ''%s -h'' to display full help.' % (PROGRAM, PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--index-file', type='str'  , dest='index_file',
                      default=os.environ['OBER'] + '/data/cgi/README.assembly_sample_subject.csv',
                      help='CGI+FINDIV ID index file')
    parser.add_option('-t', '--target-format', type='str'  , dest='format',
                      default='plink',
                      help='Target format [plink|npz]')
    parser.add_option('-c', '--recode-cgi', action='store_true', dest='recode_cgi',
                      default=False,
                      help='Recode output PLINK format as CGI letters; relevant only if -t plink is specified')
    parser.add_option('-v', '--var-file-prefix', type='str'  , dest='var_file_prefix',
                      default=os.environ['OBER'] + '/data/cgi/all.2012-09-20.testvar.chr',
                      help='tabixed CGI var file directory')
    parser.add_option('-d', '--debug'        , action='store_true'  , dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-n', '--no-swap-major-minor'        , action='store_true'  , dest='allele_swap', default=False,
                      help='Do not swap major/minor allele encoding based on the data')  # Note: reads the no-allele-swap flag and then negates it below!
    parser.add_option('-k', '--custom-identifier'        , action='store_true'  , dest='custom_id', default=False,
                      help='Assign custom identifiers to ALL snps, not just those with blank identifiers in the data')
    options, args = parser.parse_args(sys.argv[1:])
    options.allele_swap = not options.allele_swap 
    if not options.format in ['plink', 'npz']:
        print 'Unreocgnized format ''%s''. Supported formats: plink, npz' % (options.format,)
        print usage
    if len(args) != 2:
        print usage
        sys.exit(1)
    try:
        input_file = sys.stdin if args[0] == '-' else open(args[0], 'rb')
        g = extract_genotypes(input_file, index_file=options.index_file, var_file_prefix=options.var_file_prefix,
                              allele_swap=options.allele_swap, custom_id=options.custom_id, debug=options.debug)
        out = args[1]
        if options.format == 'npz':
            im.io_genotype.write(options.format, g, out + '.npz')
        elif options.format == 'plink':
            # Does not yet include the genetic map . TODO: save it to a separate PLINK file if a dedicated
            # PLINK format exists for this data
            im.io_genotype.write(options.format, g, open(out + '.tped', 'wb'),
                                 sample_id_out=out + '.tfam', recode_cgi=options.recode_cgi)
        else:
            raise ValueError('Unsupported output format ''%s''' % (options.format,))
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
