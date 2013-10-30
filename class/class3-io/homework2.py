#!/usr/bin/env python
'''
============================================================
Homework #2
- Read a PLINK map file and a corresponding frq file
into a single list with 5 columns.
- Filter the list of snp to MAF >= 5% in a bp range
[a,b]. 
- Write the result into a file.
============================================================
'''
import sys, os, csv, itertools as it
from optparse import OptionParser

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <map-file> <freq-file>\n' \
        'Read a PLINK map and freq file. Join and filter by MAF.\n' \
        'Write to standard output.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-t', '--min-maf', type='float', dest='min_maf', default=0.05,
                      help='MAF threshold - only SNPs with MAF >= this threshold are returned.')
#    parser.add_option('-l', '--bp-start', type='long', dest='bp_start', default=0,
#                      help='MAF threshold - only SNPs with MAF >= this threshold are returned.')
#    parser.add_option('-r', '--bp-stop', type='long', dest='bp_stop', default=5e8,
#                      help='MAF threshold - only SNPs with MAF >= this threshold are returned.')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 2:
        print usage
        sys.exit(1)
    
    # Read from files and filter by MAF
    with open(args[0], 'rb') as map_file:
         with open(args[1], 'rb') as frq_file:
            for map_line, frq_line in it.ifilter(lambda (map_line, frq_line): float(frq_line[4]) >= options.min_maf,
                                                 it.izip(csv.reader(map_file, delimiter='\t', skipinitialspace=True),
                                                         it.islice(csv.reader(frq_file, delimiter=' ', skipinitialspace=True), 1, None))):
                print '%2d %-30s %-10d %s/%s %f %d' % (int(map_line[0]), map_line[1], int(map_line[3]),
                                                frq_line[2], frq_line[3], float(frq_line[4]), int(frq_line[5]))
