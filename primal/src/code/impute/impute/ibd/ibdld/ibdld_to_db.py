#!/usr/bin/env python
'''
============================================================
Import IBDLD results file into the impute database.
 
Example of input line:
 1    158971    158981    22    17528721    34155457 rs5748014    
rs8141025    1567    1.663e+04    35903582 45231956    rs17750549    
rs4823813    663    9328    45909984 47347102    rs4823597    
rs5768736    213    1437    47914950 49503799    rs739193    
rs6010063    133    1589

This corresponds to the individuals (158971,158981) that share
4 IBD segments:
17528721    34155457    rs5748014    rs8141025    1567   1.663e+04
35903582    45231956    rs17750549   rs4823813    663    9328
45909984    47347102    rs4823597    rs5768736    213    1437
47914950    49503799    rs739193     rs6010063    133    1589

Created on Jul 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''

import sys, os, csv, MySQLdb, time
import numpy as np
from optparse import OptionParser

def parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    prog = os.path.basename(sys.argv[0])
    usage = 'Usage: %s [-w] [-v] [file]\n\n' \
        'Import IBDLD results from file (if one positional argument is specified) or from' \
        'standard input (if not) into the impute database.\n' \
        'Type ''%s -h'' to display full help.' % (prog, prog)

    parser = OptionParser(usage=usage)
    parser.add_option('-w', action='store_true', dest='write', default=False,
                      help='Write to database')
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    (options, args) = parser.parse_args(sys.argv[1:])

    # Argument validation
    if len(args) >= 2:
        print usage
        sys.exit(1)
        
    # Read file name
    options.file = open(args[0], 'r') if len(args) >= 1 else sys.stdin
        
    return options

def read_lines(f):
    '''Read and parse IBDLD lines from input stream f.'''
    reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
    for line in reader:
        if line:
            # Ignore tab at end of line if exists
            last_element = len(line) - 1 if line[-1] == '' else len(line)
            # Sort pair (i,j) so that i < j 
            chrom = int(line[3])
            id1 = int(min(line[1], line[2]))
            id2 = int(max(line[1], line[2]))
            # - Extract start [bp], stop [bp] for each segment
            for i in xrange(4, last_element, 6): yield chrom, id1, id2, int(line[i]), int(line[i + 1])

def print_record_counter(options, start, count):
    '''Debugging printout - #records read and elapsed time.'''
    if options.debug:
        t = time.time() - start
        print 'Read %d records in %.1f s (%.1f records/s)' % (count, t, (1.0 * count / t))

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    options = parse_command_line_args()
    
    if options.write:
        if options.debug:
            print 'Opening database connection ...\n'
        con = MySQLdb.connect(host='localhost', user='impute', passwd='impute', db='impute')
        con.autocommit(False)  # Allow bulk inserts
        cur = con.cursor()
    
    if options.debug:
        print 'Reading from %s ...' % (options.file.name,)
    # Loop over input lines, break them down to intervals and insert into our database table
    start = time.time()
    count = 0
    for data in read_lines(options.file):
        count = count + 1
        if options.debug and np.mod(count, 100000) == 0:
            print_record_counter(options, start, count)
        if options.write:
            cur.execute('insert into ibd_segment(chr,id1,id2,start_bp,stop_bp) \
            values(%d,%d,%d,%d,%d) ' % data)
    print_record_counter(options, start, count)
    
    if options.write:
        if options.debug:
            print 'Committing to database ...'
        con.commit()
        con.close()
