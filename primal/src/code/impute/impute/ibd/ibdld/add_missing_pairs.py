#!/usr/bin/env python
'''
============================================================
Add missing pairs to an IBD segment list read from standard
input. Segment list must be sorted by id1, then by id2,
with id1 <= id2.

Created on February 12, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv, optparse, numpy as np, util, traceback, itertools

#---------------------------------------------
# Constants
#---------------------------------------------
'''Dummy ID signifying EOF.'''
NO_ID = -1

#---------------------------------------------
# Methods
#---------------------------------------------
def parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    prog = os.path.basename(sys.argv[0])
    usage = 'Usage: %s [-v] <id-file>\n\n' \
        'Add missing pairs to an IBD segment list read from standard\n' \
        'input. Segment list must be sorted by id1, then by id2,\n' \
        'with id1 <= id2.\n\n' \
        'Type ''%s -h'' to display full help.' % (prog, prog)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    options, args = parser.parse_args(sys.argv[1:])

    # Argument validation
    if len(args) != 1:
        print usage
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
        
    # Read file name
    options.file = sys.stdin
    options.id_file = args[0]
    return options

def read_line(reader):
    '''Read a line from the segment input file. Return a line with dummy IDs if EOF.''' 
    try: return reader.next()
    except StopIteration: return NO_ID, NO_ID

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    try:
        options = parse_command_line_args()
        sample_id = np.sort(np.loadtxt(options.id_file).astype(np.uint))
        reader = csv.reader(options.file, delimiter=' ', skipinitialspace=True)
        line = read_line(reader)
        for i, j in ((i, j) for i, j in itertools.product(sample_id, sample_id) if i <= j):
            # Advance pointer in the list of all possible IDs to the location of the current
            # row in the file, printing empty lines till it is reached.
            if i == int(line[0]) and j == int(line[1]):
                sys.stdout.write(' '.join(line) + '\n')  # Pointer in list = current row in file 
                line = read_line(reader)  # Go to next row
            else:
                sys.stdout.write('%d %d\n' % (i, j))  # Missing pair, print empty line
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
