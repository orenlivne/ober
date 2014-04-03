#!/usr/bin/env python
'''
============================================================
Calculate call rates in imputed CGI files.

Created on February 18, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, traceback, util, cStringIO, numpy as np, impute as im
from optparse import OptionParser

MAX_HAP_TYPE = 4 # Maximum number of haplotype tag value

#---------------------------------------------
# METHODS
#---------------------------------------------
def print_count_by_snp(lines, out, variant_id, data_type, buf_size=10000, num_metadata_cols=2,
                       hap_type_low=None, hap_type_high=None):
    '''Count genotypes for each SNPs. Buffer output to the stream out. If variant_id=True, print
    variant ID as the first column of the output. If hap_type is not None, counts only genotypes
    whose hap_type tag value = hap_type.
    
    Return a list of per-sample genotype count dictionaries.'''
    # Read # samples in the first relevant line
    line = lines.next()
    num_samples, i, buf = len(line) - num_metadata_cols, 0, cStringIO.StringIO()
    # Initialize a dictionary with genotype count arrays (summed over all SNPs) of each sample
    empty_array = lambda n: [0] * n
    sample_count = dict([(x, empty_array(num_samples)) for x in im.recode.CGI_GENOTYPES])
    hap_type_filter = hap_type_low or hap_type_high
    if hap_type_filter:
        if not hap_type_low: hap_type_low = 0
        if not hap_type_high: hap_type_high = MAX_HAP_TYPE # A very big number
    while True:
        # Process line
        count = dict([(x, 0) for x in im.recode.CGI_GENOTYPES]) 
        # Count genotype values in relevant columns.
        if data_type == 'genotype':
            # Genotype data, no haplotype tagging leading character
            for sample, x in enumerate(line[num_metadata_cols:]):
                count[x] += 1
                sample_count[x][sample] += 1
        else:
            # Haplotype data, filter on haplotype type tag value if option was specified
            # Haplotypes are in the format TXX, where T=0,1,2 (phasing tag) and XX is the genotype
            if hap_type_filter:
                # If hap type does not match, lump into missing data bin ('NN')
                for sample, x in enumerate(line[num_metadata_cols:]):
                    tag = int(x[0])
                    key = x[1:] if (tag >= hap_type_low) and (tag <= hap_type_high) else 'NN'
                    count[key] += 1
                    sample_count[key][sample] += 1
            else:
                for sample, x in enumerate(line[num_metadata_cols:]):
                    key = x[1:] 
                    count[key] += 1
                    sample_count[key][sample] += 1

        # Write genotype bin counts to buffer
        if variant_id: buf.write(' '.join(line[:num_metadata_cols]) + ' ')
        for x in im.recode.CGI_GENOTYPES: buf.write('%d ' % (count[x],))
        # Write call rates: allele call rate, genotype call rate
        buf.write('%f %f\n' % (1.0 - float(count['NN']) / num_samples,
                  float(count['00'] + count['01'] + count['10'] + count['11']) / num_samples))
        
        i += 1
        if i == buf_size:
            # Write buffer to output stream, then clear buffer
            out.write(buf.getvalue())
            buf.close()
            buf = cStringIO.StringIO()
            i = 0

        # Read the next line
        try: line = lines.next()
        except StopIteration: break
    
    # Write remaining data in last buffer
    if i: out.write(buf.getvalue())
    buf.close()
    
    return sample_count

def write_sample_counts(sample_count, out):
    '''Write sample count dictionary to the file name/stream out.'''
    # Reshape data into a numpy array (rows=samples, columns=genotypes)
    num_samples, K = len(sample_count[im.recode.CGI_GENOTYPES[0]]), len(im.recode.CGI_GENOTYPES)
    count = np.zeros((num_samples, K), dtype=np.int)
    for k, x in enumerate(im.recode.CGI_GENOTYPES): count[:, k] = sample_count[x]
    np.savetxt(out, count, '%d')

####################################################################################
def __parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s\n' \
        'Calculate call rates in a CGI imputed tab-delimited standard input.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-d', '--data-file', type='str'  , dest='data_file',
                      default=None, help='If specified, reads from data file, otherwise reads from stdin')
    parser.add_option('-v', '--variant-id', action='store_true', dest='variant_id', default=False,
                      help='Output variant_id, variant_type as the first two columns')
    parser.add_option('-n', '--num-header-lines', type='int'  , dest='num_header_lines',
                      default=0, help='# header lines to skip')
    parser.add_option('-t', '--hap-type-low', type='int'  , dest='hap_type_low',
                      default=None, help='If specified, counts only genotypes whose hap_type tag value = hap_type. Values should be 1 (phased haps) or 2 (phased with known paternal origin) only. Only has effect when data_type=haplotype')
    parser.add_option('-e', '--hap-type-high', type='int'  , dest='hap_type_high',
                      default=None, help='If specified, counts only genotypes whose hap_type tag value = hap_type. Values should be 1 (phased haps) or 2 (phased with known paternal origin) only. Only has effect when data_type=haplotype')
    parser.add_option('-g', '--data-type', type='str'  , dest='data_type',
                      default='genotype', help='Data type (genotype|haplotype)')
    parser.add_option('-s', '--sample-output-file', type='str'  , dest='sample_out',
                      default='sample.out', help='Name of sample count output file.')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 0:
        print usage
        print 'No mandatory arguments should be passed in.'
        sys.exit(1)
    if options.hap_type_low and (options.hap_type_low < 0 or options.hap_type_low > MAX_HAP_TYPE):
        print usage
        sys.exit(1)
        print 'Haplotype type lower limit must be in [%d..%d].' % (0, MAX_HAP_TYPE)
    if options.hap_type_high and (options.hap_type_high < 0 or options.hap_type_high > MAX_HAP_TYPE):
        print usage
        sys.exit(1)
        print 'Haplotype type upper limit must be in [%d..%d].' % (0, MAX_HAP_TYPE)
    if not options.data_type in ['genotype', 'haplotype']:
        print usage
        print 'Unreocgnized data type ''%s''. Supported data types: genotype, haplotype' % (options.data_type,)
        sys.exit(1)

    return options

def __main(options):
    '''Main program - accepts an options struct.'''
    try:
        lines = util.read_file(open(options.data_file, 'rb') if options.data_file else sys.stdin,
                               num_header_lines=options.num_header_lines)
        # Calculate SNP counts and print to stdout; calculate sample_counts, stored in retval
        sample_count = print_count_by_snp(lines, sys.stdout, options.variant_id, options.data_type,
                                          hap_type_low=options.hap_type_low, hap_type_high=options.hap_type_high)
        # Write sample counts to sample output file
        write_sample_counts(sample_count, options.sample_out)
    except (IOError, OSError):
        traceback.print_exc(file=sys.stdout)
        sys.exit(141)

def main(**kwargs):
    '''Main program - accepts argument dictionary.'''
    # Default options
    options = util.Struct(data_file=None, variant_id=False)
    # Override with passed arguments
    options.update(**kwargs)
    # (valid, options, error_msg) = __validate_options(options)
    # if not valid:
    # raise ValueError('Bad options: %s' % (error_msg,))
    return __main(options)
    
if __name__ == '__main__':
    '''Main program - accepts CLI arguments.'''
    __main(__parse_command_line_args())
