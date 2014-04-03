#!/usr/bin/env python
'''
============================================================
Compare imputed genotypes with genotypes read from a PLINK
file (e.g., affymetrix). Calculate concordance rates. 

Created on May 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, traceback, util, impute as im, numpy as np
from optparse import OptionParser
from impute.cgi.override_imputed_by_plink import __override_imputed_by_plink

#---------------------------------------------
# Constants
#---------------------------------------------
def process_line_factory(options):
    '''A factory of process_line functors. Line processing depends on the tag value.'''
    return __process_line_concordance

def __process_line_concordance(imputed_line, g, num_metadata_cols, typed_sample_index,
                            max_flipped_hom=10, max_mismatch_partial=5, override_tag=im.constants.UNPHASED, warnings=False):
    '''Process line: override the imputed haplotypes whose indices are typed_sample_index with the
    corresponding values in the genotype file. Overridden genotypes are tagged with the tag value
    ''override_tag''. Override all genotypes g; only if h and g match and has more alleles called.'''
    H = np.array(map(list, imputed_line[num_metadata_cols:]))

    # Find sample indices that require overriding; if none found, return the original imputed line
    h = H[typed_sample_index]
    h0, h1 = h[:, 1], h[:, 2]
    g0, g1 = g[:, 0], g[:, 1]

    # Check if allele codings in h and g match; if not, flip g
    mismatch_hom = np.where(((h0 == '0') & (h1 == '0') & (g0 == 2) & (g1 == 2)) | 
                            ((h0 == '1') & (h1 == '1') & (g0 == 1) & (g1 == 1)))[0]
    flipped_hom = len(mismatch_hom)
    need_flipping = flipped_hom > max_flipped_hom
    if need_flipping: a0, a1 = '1', '0' 
    else: a0, a1 = '0', '1'

    # Total number of genotypes that are called in both data sets
    called_in_both = ((h0 != im.recode.CGI_MISSING_LETTER).astype(int) & 
                      (h1 != im.recode.CGI_MISSING_LETTER).astype(int) & 
                      (g0 != im.constants.MISSING) & 
                      (g1 != im.constants.MISSING)) 
    # Number of discordant genotypes that are called in both data sets
    discordant = called_in_both & (((h0 == a0) & (h1 == a0) & ((g0 != 1) | (g1 != 1))) | 
                                   ((((h0 == a0) & (h1 == a1)) | ((h0 == a1) & (h1 == a0))) & (g0+g1 != 3)) | 
                                   ((h0 == a1) & (h1 == a1) & ((g0 != 2) | (g1 != 2))))
    return imputed_line[:num_metadata_cols] + '%d %d %f' % (discordant, called_in_both, 1.0 - float(discordant) / called_in_both)

####################################################################################
def __parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <plink-prefix> <imputed-file>\n' \
        'Calculate concordance between imputed genotypes and a PLINK genotype data set (e.g., Affymetrix data).\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-b', '--buf-size', type='int'  , dest='buf_size',
                      default=1000, help='# lines to output at once (writing in a buffered fashion)')
    parser.add_option('-n', '--genotype-num-header-lines', type='int'  , dest='num_header_lines',
                      default=1, help='# header lines to skip in the genotype file')
    parser.add_option('-i', '--index-file', type='str'  , dest='index_file',
                      default=os.environ['OBER_DATA'] + '/cgi/README.assembly_sample_subject.csv',
                      help='CGI+FINDIV ID index file')
    parser.add_option('-m', '--num-metadata-cols', type='int'  , dest='num_metadata_cols',
                      default=8, help='# variant metadata columns preceding the genotypes in each line')
    parser.add_option('-f', '--max-flipped-hom', type='int'  , dest='max_flipped_hom',
                      default=10, help='Maximum number of allowed homozygote flips (11->22 and 22->11). If exceeds this number, plink allele coding will be flipped before overriding the imputed genotypes')
    parser.add_option('-z', '--max-partial-mismatch', type='int'  , dest='max_mismatch_partial',
                      default=5, help='Maximum number of allowed partial-genotype mismatches between original and overridden genotypes (e.g., 1N->22 or 2N->11). If exceeds this number, not overridding genptypoes at that SNP.')
    parser.add_option('-p', '--pedigree-file', type='str'  , dest='pedigree_file',
                      default=im.itu.HUTT_PED,
                      help='Pedigree file to read ALL subjects from (type+untyped, e.g., 3671 Hutterites)')
    parser.add_option('-g', '--genotype-id-file', type='str'  , dest='genotype_id_file',
                      default=im.examples.CHR22 + '/hutt.tfam',
                      help='Pedigree file to read genotyped subject IDs from (e.g., 1415 Hutterites)')
    parser.add_option('-d', '--delimiter', type='str'  , dest='delimiter', default='\\t',
                      help='Delimiter [default: %default. Use the string ''\t'' for a tab]')
    parser.add_option('-t', '--override-tag', type='int'  , dest='override_tag',
                      default=im.constants.UNPHASED, help='Tag value to tag overridden genotypes with')
    parser.add_option('-w', '--warnings', action='store_true'  , dest='warnings', default=False,
                      help='Print warning information about allele coding flipping but does not flip coding in output')
    parser.add_option('-k', '--key', type='str'  , dest='key', default='bp',
                      help='Key column to use in matching imputed and plink rows. ''bp'': \
                      match on base-pair position. ''name'': match on SNP name. Name MUST BE NUMERIC. \
                      Suitable for PLINK files generated from CGI files, where name is the CGI variant ID \
                      and uniquely identifies the SNP whereas the base-pair position does not.')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 2:
        print usage
        sys.exit(1)
    if options.delimiter == '\\t': options.delimiter = '\t'    
    return args, options

if __name__ == '__main__':
    '''Main program - accepts CLI arguments.'''
    try:
        args, options = __parse_command_line_args()
        overridden_lines = __override_imputed_by_plink(args[1], args[0], options)
        # Write overridden lines to standard output
        util.write_buffered_lines(overridden_lines, sys.stdout, buf_size=options.buf_size, delimiter=options.delimiter, process_line_factory=process_line_factory)
    except (IOError, OSError):
        traceback.print_exc(file=sys.stdout)
        sys.exit(141)
