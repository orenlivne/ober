#!/usr/bin/env python
'''
============================================================
Override imputed genotypes by affymetrix data, for those
samples and SNPs for which the latter is available. 

Created on May 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, traceback, util, impute as im, numpy as np
from optparse import OptionParser

#---------------------------------------------
# Constants
#---------------------------------------------
def key_parts(key):
    '''Return the Genotype.snp array column name and the corresponding
    imputed file field number for the key type ''key''.'''
    if key == 'bp': return 'base_pair', 3
    elif key == 'name': return 'name', 0 
    else: raise ValueError('Unsupported key type ''%s''' % (key,))

def process_lines(process_line, imputed_lines, g, num_metadata_cols, typed_sample_index,
                  max_flipped_hom=10, max_mismatch_partial=5,
                  override_tag=im.constants.UNPHASED, warnings=False, key='bp'):
    '''Override genotypes and yield the line of new genotype items for each SNP.'''
    # List of base-pairs of the SNPs in the Genotype object g, mapped to the index in the g.snp array (and therefore also in the g.data array)
    snp_field_name, field_num = key_parts(key)
    f = g.snp[snp_field_name]
    if len(f) == 1: f = [str(f[0])]
    typed_snp_index = dict((int(v), k) for k, v in enumerate(f))
    for line in imputed_lines:
        # Process line
        key = int(line[field_num])
        if key in typed_snp_index:
            # Imputed SNP found in PLINK file, override genotypes when possible
            overridden_line = process_line(line, g.data_of_snp(typed_snp_index[key]), num_metadata_cols,
                                           typed_sample_index, override_tag=override_tag, warnings=warnings,
                                           max_flipped_hom=max_flipped_hom, max_mismatch_partial=max_mismatch_partial)
        else: overridden_line = line
        yield overridden_line

def process_line_factory(options):
    '''A factory of process_line functors. Line processing depends on the tag value.'''
    if options.override_tag == im.constants.LD_IMPUTED_WITH_ORIGIN:
        return __process_line_ld_poo
    else:
        return __process_line_priority

def __process_line_priority(imputed_line, g, num_metadata_cols, typed_sample_index,
                            max_flipped_hom=10, max_mismatch_partial=5, override_tag=im.constants.UNPHASED, warnings=False):
    '''Process line: override the imputed haplotypes whose indices are typed_sample_index with the
    corresponding values in the genotype file. Overridden genotypes are tagged with the tag value
    ''override_tag''. Override all genotypes g; only if h and g match and has more alleles called.'''
    H = np.array(map(list, imputed_line[num_metadata_cols:]))

    # Find sample indices that require overriding; if none found, return the original imputed line
    h = H[typed_sample_index]
    h0, h1 = h[:, 1], h[:, 2]
    g0, g1 = g[:, 0], g[:, 1]
    need_override = np.where((h0 != im.recode.CGI_MISSING_LETTER).astype(int) + (h1 != im.recode.CGI_MISSING_LETTER).astype(int) < 
                             (g0 != im.constants.MISSING).astype(int) + (g1 != im.constants.MISSING).astype(int))[0]
    if not need_override.size: return imputed_line
        
    # Need to override; check if allele codings in h and g match; if not, flip g
    mismatch_hom = np.where(((h0 == '0') & (h1 == '0') & (g0 == 2) & (g1 == 2)) | 
                             ((h0 == '1') & (h1 == '1') & (g0 == 1) & (g1 == 1)))[0]
    flipped_hom = len(mismatch_hom)
    need_flipping = flipped_hom > max_flipped_hom
    if warnings and need_flipping:
        # If warning mode, don't override if allele coding needs to be flipped
        sys.stderr.write('=== WARNING: allele coding mismatch, not overriding. Metadata %s flipped_hom %d mismatch_hom %s\n' % ('\t'.join(imputed_line[:num_metadata_cols]), flipped_hom, repr(mismatch_hom)))
        return imputed_line
    if need_flipping: a0, a1 = '1', '0' 
    else: a0, a1 = '0', '1'

    # Check for partial genotype mismatches: h=22 and g=1N, or h=11 and g=0N
    mismatch_partial = np.where((((h0 == a0) | (h1 == a0)) & (g0 == 2) & (g1 == 2)) | 
                                (((h0 == a1) | (h1 == a1)) & (g0 == 1) & (g1 == 1)))[0]
    if len(mismatch_partial) > max_mismatch_partial:
        if warnings: sys.stderr.write('=== WARNING: too many partial genotype mismatches found %d, not overriding. Mismatches %s\n' % (len(mismatch_partial), repr(mismatch_partial)))
        return imputed_line
        
    t = typed_sample_index[need_override]
    H[t, 1:] = im.recode.recode_cgi_flipped(g[need_override]) if need_flipping else im.recode.recode_cgi(g[need_override])  
    H[t, 0] = override_tag
    return imputed_line[:num_metadata_cols] + [str(h[0]) + str(h[1]) + str(h[2]) for h in H]

def __process_line_ld_poo(imputed_line, g, num_metadata_cols, typed_sample_index,
                           max_flipped_hom=10, max_mismatch_partial=5, override_tag=im.constants.UNPHASED, warnings=False):
    '''Process line: override LD-based genotypes h by LD-based PO-haplotypes g (note: the h, g variable
    roles are reversed w.r.t. __process_line_ld_poo()). Override only if both
    are called and they match (except possibly for allele ordering).'''
    H = np.array(map(list, imputed_line[num_metadata_cols:]))
    # Find sample indices that require overriding; if none found, return the original imputed line
    h = H[typed_sample_index]
    old_tag, h0, h1 = h[:, 0].astype(int), h[:, 1], h[:, 2]
    # Prevent error converting 'N' to int
    h0[h0 == im.recode.CGI_MISSING_LETTER] = '9'  # Large dummy value
    h1[h1 == im.recode.CGI_MISSING_LETTER] = '9'  # Large dummy value
    np.set_printoptions(threshold=np.nan)
    g0, g1 = g[:, 0], g[:, 1]
    # Override only 3xx genotypes; override only if both g, h are fully called and their dosages match
    # Note: h alleles are 0-based, g are 1-based (need to account for that when comparing dosages) 
    need_override = np.where((old_tag == im.constants.LD_IMPUTED)
                             & (h0 != im.recode.CGI_MISSING_LETTER) & (h1 != im.recode.CGI_MISSING_LETTER)
                             & (g0 != im.constants.MISSING) & (g1 != im.constants.MISSING)
                             & (h0.astype(int) + h1.astype(int) + 2 == g0 + g1))[0]
    if not need_override.size: return imputed_line
    # Override h (genotypes) by g (haplotypes)
    t = typed_sample_index[need_override]
    H[t, 1:] = im.recode.recode_cgi(g[need_override])  
    H[t, 0] = override_tag
    return imputed_line[:num_metadata_cols] + [str(h[0]) + str(h[1]) + str(h[2]) for h in H]

def __override_imputed_by_plink(imputed_file, plink_prefix, options, process_line_factory=process_line_factory):
    '''Main program for overriding imputing genotype file by a PLINK data file. Accepts an options struct.'''
    # Read plink data
    g = im.io_genotype.read('plink', 'genotype', prefix=plink_prefix, lazy_load=True)
    # Read WGS sample indices from the pedigree. Must match the genotype files' ordering.
    typed_sample_index = np.array(map(im.io_pedigree.read(options.pedigree_file, options.genotype_id_file)._node_of.get, g.sample_id))
    # Read imputed line pairs         
    imputed_lines = util.read_file(imputed_file, num_header_lines=options.num_header_lines, delimiter=options.delimiter)
    # Override lines by data from g (depending on the scenario specified by options)
    process_line = process_line_factory(options)
    return process_lines(process_line, imputed_lines, g, options.num_metadata_cols, typed_sample_index,
                         max_flipped_hom=options.max_flipped_hom, max_mismatch_partial=options.max_mismatch_partial,
                         override_tag=options.override_tag, warnings=options.warnings, key=options.key)

def override_imputed_by_plink(imputed_file, plink_prefix, process_line_factory=process_line_factory, **kwargs):
    '''Main program for overriding imputing genotype file by a PLINK data file. Accepts an argument dictionary.'''
    # Default options
    options = util.Struct(num_header_lines=0, num_metadata_cols=8, override_tag=im.constants.UNPHASED,
                          max_flipped_hom=10, max_mismatch_partial=5, delimiter='\t',
                          index_file=os.environ['OBER_DATA'] + '/cgi/README.assembly_sample_subject.csv',
                          pedigree_file=im.itu.HUTT_PED, genotype_id_file=im.examples.CHR22 + '/hutt.tfam',
                          warnings=False, key='bp')
    # Override with passed arguments
    options.update(**kwargs)
    return __override_imputed_by_plink(imputed_file, plink_prefix, options, process_line_factory=process_line_factory)
    
####################################################################################
def __parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <plink-prefix> <imputed-file>\n' \
        'Override CGI imputed tab-delimited at a set of SNPs with the genotypes read a plink file.\n' \
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
        util.write_buffered_lines(overridden_lines, sys.stdout, buf_size=options.buf_size, delimiter=options.delimiter)
    except (IOError, OSError):
        traceback.print_exc(file=sys.stdout)
        sys.exit(141)
