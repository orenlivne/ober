#!/usr/bin/env python
import os, sys, csv, itertools as it
from optparse import OptionParser

def parse_command_line(argv):
    PROGRAM = os.path.basename(argv[0])
    usage = 'Usage: %s <map-file> <freq-file> <out-file>\n' \
        'Read a PLINK map and freq file. Join and filter by MAF.\n' \
        'Write to standard output.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-t', '--min-maf', type='float', dest='min_maf', default=0.05,
                      help='MAF threshold - only SNPs with MAF >= this threshold are returned.')
    options, args = parser.parse_args(argv[1:])
    if len(args) != 3:
        print usage
        sys.exit(1)
    return options, args

def process_files(map_file_name, frq_file_name, out_file_name, min_maf=0.05):
    try:
        with open(map_file_name, 'rb') as map_file:
            with open(frq_file_name, 'rb') as frq_file:
                with open(out_file_name, 'wb') as out_file:
                    for map_line, frq_line in it.ifilter(lambda (_, y): float(y[4]) >= min_maf,
                                                        zip(
                                                            csv.reader(map_file, delimiter='\t', skipinitialspace=True),
                                                            it.islice(csv.reader(frq_file, delimiter=' ', skipinitialspace=True), 1, None)
                                                            )):
                        out_file.write('%2d %-20s %-10d %s/%s %.3f\n' % \
                                       (int(map_line[0]), map_line[1], int(map_line[3]), frq_line[2], frq_line[3], float(frq_line[4])))
    except IOError as e:
        print e
    
if __name__ == '__main__':
    options, args = parse_command_line(sys.argv)
    process_files(args[0], args[1], args[2], min_maf=options.min_maf)
