#!/usr/bin/env python
# encoding: utf-8
# filename: count_cgi_genotypes.py
# cython: profile=True
'''
============================================================
Calculate call rates in imputed CGI files.

Created on February 18, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv, itertools, numpy as np, traceback, util
from optparse import OptionParser

#---------------------------------------------
# Constants
#---------------------------------------------
# All possible genotypes
GENOTYPES = [x[0] + x[1] for x in list(itertools.product('N01', 'N01'))]
# Converts CGI allele code to our numeric code
CGI_LETTER_TO_ALLELE = {'N': 0, '0': 1, '1': 2}

def genotype_start_index(line):
    '''Return the start index of g entries in the list line. If not found, returns -1.'''
    index = 6
    for x in line[6:]:
        if x in GENOTYPES:
            return index
        else:
            index += 1
    return -1

def print_count_by_snp(lines, out, id_list):
    '''Count total genotypes for each SNPs.'''
    # Initialize all genotype counts at 0
    # Stream lines and increment counts
    for line in lines:
        # Lines may start with a variable no. of items from the csv reader's perspective (e.g.,
        # indel with empty substitution fields will result in consecutive spaces. Calculate the
        # start of the genotype sublist  
        index = genotype_start_index(line)
        genotype = line[index:]
        # Pick out the relevant IDs
        count = dict(zip(GENOTYPES, [0] * len(GENOTYPES)))
        for x in (np.array(genotype)[id_list] if id_list is not None else genotype):
            count[x] += 1
        print_count_total(count, out)

def count_total(lines, id_list):
    '''Count total genotypes over the entire file.'''
    # Initialize all genotype counts at 0
    count = dict(zip(GENOTYPES, [0] * len(GENOTYPES)))
    print count
    # Stream lines and increment counts
    for line in lines:
        # Lines may start with a variable no. of items from the csv reader's perspective (e.g.,
        # indel with empty substitution fields will result in consecutive spaces. Calculate the
        # start of the genotype sublist  
        index = genotype_start_index(line)
        genotype = line[index:]
        # Pick out the relevant IDs
        for x in (np.array(genotype)[id_list] if id_list is not None else genotype):
            count[x] += 1
    return count
    
def print_count_total(count, out):
    '''Print total count results: (genotype count frequency) columns for all genotypes.'''
    total = sum(count.itervalues())
    for k in GENOTYPES:
        out.write('%s %8d %.3f ' % (''.join(map(str, map(CGI_LETTER_TO_ALLELE.get, k))), count[k], (1.0 * count[k]) / total))
    out.write('\n')

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
    parser.add_option('-i', '--id-index-file', type='str'  , dest='id_file',
                      default=None, help='If specified, outputs only the IDs listed in this file (these are indices between 0 and #ids-1, if the input file has #ids genotype columns)')
    parser.add_option('-s', '--snp', action='store_true'  , dest='group_by_snp', default=False,
                      help='Group by snp')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 0:
        print usage
        sys.exit(1)
    return options

def __main(options):
    '''Main program - accepts an options struct.'''
    
    # If id file is specified, read into the 'id_list' array
    id_list = np.loadtxt(options.id_file, dtype=np.int) if options.id_file else None 
    
    # Init dictionary of all-possible-genotypes-to-counts
    try:
        f = open(options.data_file, 'rb') if options.data_file else sys.stdin
        lines = (line for line in csv.reader(f, delimiter='\t', skipinitialspace=True) if line)
        if options.group_by_snp:
            print_count_by_snp(lines, sys.stdout, id_list)
        else:
            count = count_total(lines, id_list)
            print_count_total(count, sys.stdout)
    except (IOError, OSError):
        traceback.print_exc(file=sys.stdout)
        sys.exit(141)

def main(**kwargs):
    '''Main program - accepts argument dictionary.'''
    # Default options
    options = util.Struct(data_file=None, id_file=None, group_by_snp=False)
    # Override with passed arguments
    options.update(**kwargs)
    # (valid, options, error_msg) = __validate_options(options)
    # if not valid:
    # raise ValueError('Bad options: %s' % (error_msg,))
    return __main(options)
    
if __name__ == '__main__':
    '''Main program - accepts CLI arguments.'''
    __main(__parse_command_line_args())
