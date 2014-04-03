#!/usr/bin/env python
'''
============================================================
Override imputed genotypes of the 98 WGS Hutterites by
genotype data, if it is available.

Created on February 18, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv, itertools, traceback, util, cStringIO, impute as im, db_gene, numpy as np
from optparse import OptionParser

#---------------------------------------------
# Constants
#---------------------------------------------
# All possible genotypes
MISSING = 'N'
GENOTYPES = [x[0] + x[1] for x in list(itertools.product(MISSING + '01', MISSING + '01'))]

def process_lines(genotype_lines, imputed_lines, out, num_metadata_cols, wgs_sample_index,
                  buf_size=10000, delimiter='\t'):
    '''Count genotypes for each SNPs. Buffer output to the stream out. If variant_id=True, print
    variant ID as the first column of the output.'''
    # Read # samples in the first relevant line
    genotype_line = genotype_lines.next()
    imputed_line = imputed_lines.next()
    i, buf = 0, cStringIO.StringIO()
    
    while True:
        i += 1
        if i == buf_size:
            # Write buffer to output stream, then clear buffer
            out.write(buf.getvalue())
            buf.close()
            buf = cStringIO.StringIO()
            i = 0
        process_line(genotype_line, imputed_line, num_metadata_cols, wgs_sample_index, buf)
        # Read the next line
        try:
            genotype_line = genotype_lines.next()
            imputed_line = imputed_lines.next()
        except StopIteration:
            break
    if i:
        out.write(buf.getvalue())
    buf.close()

def process_line(genotype_line, imputed_line, num_metadata_cols, wgs_sample_index, buf, delimiter='\t'):
    '''Process line: overrwite the imputed genotypes whose indices are wgs_sample_index with the
    corresponding values in the genotype file.'''
    buf.write(delimiter.join(genotype_line[:num_metadata_cols]))
    G = np.array(map(list, genotype_line[num_metadata_cols:]))
    H = np.array(map(list, imputed_line[num_metadata_cols:]))
    # Override imputed by genotype, if necessary
    h = H[wgs_sample_index]
    need_override = np.where((h[:, 0] != 'N').astype(int) + (h[:, 1] != 'N') < 
                             (G[:, 0] != 'N').astype(int) + (G[:, 1] != 'N'))[0]
    if need_override.size:
        # print H[wgs_sample_index[need_override]]
        H[wgs_sample_index[need_override]] = G[need_override]
        # print H[wgs_sample_index[need_override]]
        
    # Write updated imputed genotypes to buffer
    buf.write(delimiter + delimiter.join(str(h[0]) + str(h[1]) for h in H) + '\n')
        
def read_file(file_name, num_header_lines=0, delimiter='\t'):
    '''Return a generator expression of the lines in the file file_name, skipping the first num_header_lines.'''
    lines = (line for line in csv.reader(open(file_name, 'rb'), delimiter=delimiter, skipinitialspace=True) if line)
    # Skip header lines
    for _ in xrange(num_header_lines): lines.next()
    return lines 

####################################################################################
def __parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <genotype-file> <imputed-file>\n' \
        'Calculate call rates in a CGI imputed tab-delimited standard input.\n' \
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
    parser.add_option('-p', '--pedigree-file', type='str'  , dest='pedigree_file',
                      default=im.itu.HUTT_PED,
                      help='Pedigree file to read ALL subjects from (type+untyped, e.g., 3671 Hutterites)')
    parser.add_option('-g', '--genotype-id-file', type='str'  , dest='genotype_id_file',
                      default=im.examples.CHR22 + '/hutt.tfam',
                      help='Pedigree file to read genotyped subject IDs from (e.g., 1415 Hutterites)')
    parser.add_option('-d', '--delimiter', type='str'  , dest='delimiter', default='\\t',
                      help='Delimiter [default: %default. Use the string ''\t'' for a tab]')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 2:
        print usage
        sys.exit(1)
    options.genotype_file = args[0]
    options.imputed_file = args[1]
    if options.delimiter == '\\t': options.delimiter = '\t'    
    return options

def __main(options):
    '''Main program - accepts an options struct.'''
    try:
        # Read WGS sample indices from the pedigree. Must match the genotype files' ordering.
        wgs_sample_index = np.array(map(im.io_pedigree.read(options.pedigree_file, options.genotype_id_file)._node_of.get, db_gene.cgi.ids.cgi_id_to_sample_id(options.index_file).values()))

        # Read and process (genotype, imputed) line pairs         
        genotype_lines = read_file(options.genotype_file, num_header_lines=options.num_header_lines, delimiter=options.delimiter)
        imputed_lines = read_file(options.imputed_file, delimiter=options.delimiter)
        process_lines(genotype_lines, imputed_lines, sys.stdout, options.num_metadata_cols, wgs_sample_index, delimiter=options.delimiter, buf_size=options.buf_size)
    except (IOError, OSError):
        traceback.print_exc(file=sys.stdout)
        sys.exit(141)

def main(**kwargs):
    '''Main program - accepts argument dictionary.'''
    # Default options
    options = util.Struct(num_header_lines=0, num_metadata_cols=8, buf_size=10000,
                          pedigree_file=im.itu.HUTT_PED, genotype_file=im.examples.CHR22)
    # Override with passed arguments
    options.update(**kwargs)
    # (valid, options, error_msg) = __validate_options(options)
    # if not valid:
    # raise ValueError('Bad options: %s' % (error_msg,))
    return __main(options)
    
if __name__ == '__main__':
    '''Main program - accepts CLI arguments.'''
    __main(__parse_command_line_args())
