#!/usr/bin/env python
'''
============================================================
Extract imputed Hutterite genotypes from a tabixed imputed
genotype file set.

Prerequisites: tabix

Created on April 17, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv, util, StringIO, traceback, re, itertools as it, numpy as np
from optparse import OptionParser
from impute.impute_test_util import HUTT_PED
from impute.phasing.examples import CHR22
from impute.tools.recode import recode12
from impute.data import io_pedigree
from db_gene.ucsc.ucsc_dao import SnpDao

#---------------------------------------------
# Constants
#---------------------------------------------
# Path to tabix program
TABIX_CMD = 'tabix' 

#---------------------------------------------
# Methods
#---------------------------------------------
def __parse_command_line_args(argv):
    '''Parse and validate command-line arguments.'''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(argv[0])
    usage = 'Usage: %s [options] <var_file_prefix> <input_file> <output_file>\n' \
        'Extract CGI or imputed genotypes from tabixed genotype files var_file_prefix.chr*.gz for a\n' \
        'list of locations/SNP RS numbers. Reads from stdin. Each stdin line should be\n' \
        '\n either an rs number, or a base-pair location or segment in the format chrX:YYY[-ZZZ],\n' \
        'where X=chromosome, YYY=base-pair start, ZZZ=base-pair end, not inclusive\n' \
        '(if not specified, ZZZ=YYY+1).\n\n' \
        'Supports autosomal chromosomes only.\n' \
        '\n' \
        'Type ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-p', '--pedigree-file', type='str'  , dest='pedigree_file',
                      default=HUTT_PED,
                      help='Pedigree file to read ALL subjects from (type+untyped, e.g., 3671 Hutterites)')
    parser.add_option('-g', '--genotype-id-file', type='str'  , dest='genotype_id_file',
                      default=CHR22 + '/hutt.tfam',
                      help='Pedigree file to read genotyped subject IDs from (e.g., 1415 Hutterites)')
    parser.add_option('-d', '--db-url', type='str'  , dest='db_url',
                      default='mysql://ucsc:ucsc@localhost/ucsc',
                      help='URL of UCSC database for locating RS numbers, if ''rs'' format is specified.')
    parser.add_option('-i', '--id-index-file', type='str'  , dest='id_file',
                      default=None, help='If specified, outputs only the FINDIVs listed in this file (one per line)')
    parser.add_option('-f', '--genotype-filter', type='str'  , dest='genotype_filter',
                      default='all',
                      help='Genotype filter to apply, if output format = ''lgen''. all=no filter (all genotypes included). full=fully-called genotypes only.')
    parser.add_option('-o', '--output-format', type='str'  , dest='output_format',
                      default='cgi',
                      help='Output format. cgi=3-letter imputed genotypes in (extended) CGI format. matrix=one line per snp with all genotypes. lgen=LGEN format: one file for all variants whose columns are snp_id sample_id allele1 allele2. tped=PLINK transposed format')
    parser.add_option('-l', '--genotype-format', type='str'  , dest='genotype_format',
                      default='num',
                      help='Genotype Output format. num=numerical (allele coding: 0/1/N). letter=ref allele/minor allele/missing.)')
    parser.add_option('-k', '--custom-identifier'        , action='store_true'  , dest='custom_id', default=False,
                      help='Assign custom identifiers to ALL snps, not just those with blank identifiers in the data')
    parser.add_option('-v', '--debug'        , action='store_true'  , dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-r', '--print-frequency'        , type='int'  , dest='print_frequency',
                      default=100, help='Print a progress message every these many records. If non-positive, progress printouts are suppressed. If negative, all printouts are suppressed.')
    parser.add_option('-b', '--break-ties-to-snps'        , action='store_true'  , dest='break_ties_to_snps', default=False,
                      help='If multiple matches are found, restrict them to SNPs only')
    parser.add_option('-s', '--single-matches'        , action='store_true'  , dest='single_matches', default=False,
                      help='Exclude variants with multiple matches in the CGI data')
    parser.add_option('-n', '--no-na-calls'        , action='store_true'  , dest='remove_na', default=False,
                      help='Exclude nameless variants (i.e. whose name in CGI is ''-'')')
    parser.add_option('-t', '--type', type='str'  , dest='type',
                      default='imputed', help='data type (genotype=CGI genotype - 2 alleles; imputed=imputed - 1 tag + 2 alleles)')
    parser.add_option('-y', '--variant-type', type='str'  , dest='variant_type',
                      default=None,
                      help='Output only variants of this type (snp|ins|del|sub). If not specified, all variants are output.')
    options, args = parser.parse_args(argv[1:])
    if not options.type in ['genotype', 'imputed']:
        print 'Unreocgnized input format ''%s''. Supported formats: genotype, imputed' % (options.output_format,)
        print usage
        sys.exit(1)
    if not options.output_format in ['cgi', 'matrix', 'lgen', 'tped']:
        print 'Unreocgnized output format ''%s''. Supported formats: cgi, matrix, lgen' % (options.output_format,)
        print usage
        sys.exit(1)
    if not options.genotype_format in ['num', 'letter']:
        print 'Unreocgnized genotype output format ''%s''. Supported formats: num, letter' % (options.genotype_format,)
        print usage
        sys.exit(1)
    if not options.genotype_filter in ['all', 'full']:
        print 'Unreocgnized variant format ''%s''. Supported formats: all, full' % (options.genotype_filter,)
        print usage
        sys.exit(1)
    if options.variant_type and not options.variant_type in ['snp', 'ins', 'del', 'sub']:
        print 'Unreocgnized variant type ''%s''. Supported formats: snp, ins, del, sub' % (options.variant_type,)
        print usage
        sys.exit(1)
    if len(args) != 3:
        print usage
        sys.exit(1)
    return options, args

__SNP_CHR_REG = re.compile('^chr(\d+)$')
__CHR_REGEX = re.compile('chr(\d+):(\d+)(.*)')
__RS_REGEX = re.compile('(rs.*)')
__NOT_FOUND = 0

def __parse_line(line, snp_dao, debug=False):
    '''Parse a line into a list of locations tuple (chrom, (bp_start, bp_stop)).''' 
    m = __CHR_REGEX.match(line)
    if m:
        # Location specified in the format chrX:YYY[-ZZZ]
        chrom, bp_start = int(m.group(1)), int(m.group(2))
        bp_stop = int(m.group(3)[1:]) if m.group(3) else (bp_start + 1)
        retval = [(line, chrom, (bp_start, bp_stop))]
        if debug: print line, retval
        return retval
    
    if snp_dao.online:
        m = __RS_REGEX.match(line)
        if m:
            # RS # specified. Exclude names like 'chr6_apd_hap1'.
            matches = [(line, int(k.group(1)), (snp.bp, snp.bp))
                       for (snp, k) in ((snp, __SNP_CHR_REG.match(snp.chrom)) for snp in snp_dao.get_snps([line])) if k]
            retval = matches if matches else [(line, __NOT_FOUND, (line, None))]
            if debug: print retval
            return retval
    
    raise ValueError('Unrecognized location format: ''%s''' % (line,))

'''Converts a genotype filter string into a filter functor.'''
__GENOTYPE_FILTER = {'full': lambda g: g[0] != 'N' and g[1] != 'N',
                     'all': lambda g: True}

'''Removes partial genotype calls, if filter is specified. m=missing allele code'''
__GENOTYPE_CLEANER = {'full': lambda g, m: m + m if ((g[0] != m) ^ (g[1] != m)) else g,
                     'all': lambda g, m: g}

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    options, args = __parse_command_line_args(sys.argv)
    var_file_prefix, input_file, output_file = args
    genotype_filter = __GENOTYPE_FILTER[options.genotype_filter]
    genotype_cleaner = __GENOTYPE_CLEANER[options.genotype_filter]
    allele_start_index = 1 if options.type == 'imputed' else 0
    try:
        # Read location list
        snp_dao = SnpDao(options.db_url)
        input_file = sys.stdin if args[1] == '-' else open(input_file, 'rb')
        
        # Sample IDs are read from the pedigree, and must match the genotype files' ordering
        pedigree = io_pedigree.read(options.pedigree_file, options.genotype_id_file)
        node_of = pedigree.node_of
        sample_id = np.loadtxt(options.id_file, dtype=np.int, ndmin=1) if options.id_file else pedigree._sample_id
        N = pedigree.num_genotyped
        sample_index = np.array(filter(lambda x: x is not None and x < N, map(node_of.get, sample_id)))  # Filter FINDIVs that are not in the imputed FINDIV set
                          
        out = open(output_file, 'wb')
        # Print header line
        if options.output_format == 'matrix':
            out.write('\t'.join(['variant', 'chromosome', 'bp_start', 'bp_stop', 'variant_type', 'ref_allele', 'minor_allele'] + map(str, sample_id)) + '\n')
                
        # Extract data using tabix for each location. A location may be a range of bps and may
        # correspond to multiple output lines.
        snp_count = [0, 0, 0, 0, 0]  # #found variants; #not-found variants; #multiply-matching variants; #nameless variants
#        print list(enumerate(it.chain.from_iterable(__parse_line(line, snp_dao, debug=options.debug) 
#                                                                                                for line in (line.rstrip('\n').rstrip('\r') 
#                                                                                                             for line in input_file) if line), 1))
        for row, (input_line, chrom, (bp_start, bp_stop)) in enumerate(it.chain.from_iterable(__parse_line(line, snp_dao, debug=options.debug) 
                                                                                                for line in (line.rstrip('\n').rstrip('\r') 
                                                                                                             for line in input_file) if line), 1):
            if options.debug:
                print '-' * 70
                print 'input_line', input_line
            if options.print_frequency > 0 and row % options.print_frequency == 0: print 'Processed %d input lines' % (row,)
            if chrom == __NOT_FOUND:
                # No match for RS# / location in SNP database
                found = False
                print 'Did not find a match for ''%s'' in SNP database' % (input_line,)
            else:
                # No match for RS# / location in genotype file
                cmd = '%s %s.chr%d.tsv.gz chr%d:%d-%d' % (TABIX_CMD, var_file_prefix, chrom, chrom, bp_start, bp_stop)
                if options.debug:
                    print 'Running command', cmd
                _, output, _ = util.run_command(cmd, verbose=options.debug)
                if options.debug:
                    print 'output', output 
                    name = bp_start
                lines = list(csv.reader(StringIO.StringIO(output), skipinitialspace=True, delimiter='\t'))
                    
                # If there are multiple matches for a point-location, exclude lines that map to multiple
                # dbsnp SNP identifiers (separated by ;)
                if len(lines) > 1:
                    if options.variant_type: lines = filter(lambda x: line[4] != options.variant_type, lines)
                    if options.single_matches: lines = filter(lambda x: not ';' in x[7], lines)
                    if options.remove_na and not options.custom_id: lines = filter(lambda x: x[7] and x[7] != '-', lines)
                    # If multiple lines and the tie-breaking flag is specified, only return the SNP results
                    if options.break_ties_to_snps: lines = filter(lambda x: x[4] == 'snp', lines)
                    
                found = len(lines) > 0
                if not found: print 'Did not find a match for ''%s'' (chr%d:%d-%d) in CGI files' % (input_line, chrom, bp_start, bp_stop)
            
            if not found:
                # Variant not found, output a dummy line
                if options.output_format == 'matrix': out.write(name + '\t' + '\t'.join(['N/A'] * (len(sample_id) + 7)) + '\n')
                else: pass  # LGEN format: don't output anything
            else:
                # Variant found, output real data
                if len(lines) > 1:
                    print 'Multiple (%d) matches for ''%s'' (chr%d:%d-%d) in CGI files' % (len(lines), input_line, chrom, bp_start, bp_stop)
                    snp_count[2] += 1
                    if options.single_matches: continue
                for line in lines:
                    # If SNP name is missing, create a unique identifier by concatenating the chromosome, end position and minor allele letter
                    name = line[7] if line[7] and not options.custom_id else 'chr%d_%d_%s' % (chrom, int(line[3]), line[6])
                    if name == '-':
                        # Filter nameless variants
                        if options.remove_na:
                            snp_count[3] += 1
                            continue 
                        name = 'N/A'

                    # Filter on variant type
                    if options.variant_type and (line[4] != options.variant_type):
                        print 'Wrong variant type ''%s'' (chr%d:%d-%d) in CGI files' % (line[4], chrom, bp_start, bp_stop)
                        snp_count[4] += 1
                        continue

                    # Autosomes - imputed genotype data (tag + two letters per sample. Read the letters.)
                    genotypes = np.array([[item[allele_start_index], item[allele_start_index + 1]] for item in line[8:]])
                    sample_id_filtered = sample_id.copy()
                    if options.debug:
                        print 'genotypes', genotypes
                        print 'sample_id_filtered', sample_id_filtered
                        
                    # Filter samples
                    if options.id_file:
                        genotypes = genotypes[sample_index]
                    if options.output_format == 'lgen' and options.genotype_filter != 'all':
                        filtered = np.where(np.array(map(genotype_filter, genotypes)))[0]
                        genotypes = genotypes[filtered]
                        sample_id_filtered = sample_id_filtered[filtered]
                    
                    # Convert genotypes to desired format
                    if options.genotype_format == 'letter':
                        allele_letter = {'0': line[5], '1': line[6], 'N': '0' if options.output_format == 'tped' else '.'}
                        if options.debug:
                            print 'allele_letter', allele_letter
                        genotypes = [[allele_letter[g[0]], allele_letter[g[1]]] for g in genotypes]
                    
                    # Write genotypes to a separate file or to the same file depending on the output format
                    if options.output_format == 'cgi':
                        # Same as input format, just output the requested subset of the samples 
                        out.write('\t'.join(line[:8]) + '\t' +
                                  '\t'.join(line[i + 8] for i in sample_index) + '\n')
                    elif options.output_format == 'matrix':
                        # Matrix of genotypes (SNP names are not printed)
                        out.write(name + '\t' + '\t'.join(line[1:7]) + '\t' + 
                                  '\t'.join(g[0] + '\\' + g[1] for g in genotypes) + '\n')
                    elif options.output_format == 'lgen':
                        # Listing-by-long format (PLINK LGEN)
                        if options.debug: sys.stdout.write('\n'.join('%s\t%d\t%s\t%s' % (name, sample_id_filtered[i], g[0], g[1]) for i, g in enumerate(genotypes)))
                        out.write('\n'.join('%s\t%d\t%s\t%s' % (name, sample_id_filtered[i], g[0], g[1]) for i, g in enumerate(genotypes)) + '\n')
                    elif options.output_format == 'tped':
                        # PLINK transposed format
                        if options.genotype_format == 'num': recoder = recode12
                        else: recoder = lambda g: str(g[0]) + ' ' + str(g[1])
                        missing_code = '0' if options.genotype_format == 'letter' else 'N'
                        out.write('%d %s %d %d ' % (chrom, name, 0, bp_start) + 
                                  ' '.join(map(recoder, map(lambda g: genotype_cleaner(g, missing_code), genotypes))) + '\n')

            snp_count[0 if found else 1] += 1
        if options.output_format == 'matrix': 
            out.close()
        if options.print_frequency >= 0:
            print 'Requested %d locations: %d found, %d not found, %d duplicates, %d nameless, %d other types.' % (row, snp_count[0], snp_count[1], snp_count[2], snp_count[3], snp_count[4])
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
