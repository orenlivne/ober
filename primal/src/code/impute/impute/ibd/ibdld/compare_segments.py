#!/usr/bin/env python
'''
============================================================
Compare two segment files (e.g., from our program vs.
IBDLD) in standard format. Output segment set concordance
measures for each sample pair.

Created on February 13, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv, optparse, impute as im, util, traceback

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
    usage = 'Usage: %s [-v] <file1> <file2>\n\n' \
    'Compare two segment files (e.g., from our program vs.\n' \
    'IBDLD) in standard format. Output segment set concordance\n' \
    'measures for each sample pair.\n\n' \
    'Type ''%s -h'' to display full help.' % (prog, prog)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    options, args = parser.parse_args(sys.argv[1:])
    
    # Argument validation
    if len(args) != 2:
        print usage
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
        
    # Read file name
    options.file1 = args[0]
    options.file2 = args[1]
    return options

def parse_line(line):
    '''Parse a segment file line into a key (sample1,sample2) and value (segment set).'''
    items = [int(x) for x in line]
    return tuple(items[0:2]), im.segment.DisjointSegmentSet(zip(items[2::2], items[3::2])) 

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    options = parse_command_line_args()
    file_reader = lambda f: csv.reader(open(f, 'rb'), delimiter=' ', skipinitialspace=True)
    param = im.PhaseParam()
    try:
        reader1 = file_reader(options.file1)
        reader2 = file_reader(options.file2)
        for i, line in enumerate(reader1):
            key, A = parse_line(line)
            key2, B = parse_line(reader2.next())
            if key != key2:
                raise ValueError('Sample pairs are not the same in both files: (%d,%d), (%d,%d) at line %d' % \
                                 key + key2 + (i + 1,))
            else:
                # Output statistics. Only log pairs for which samples were found
                len_A, len_B = A.length, B.length
                if len_A != 0 or len_B != 0:
                    f = param.kinship(key[0], key[1])
                    sys.stdout.write('%d %d %f %d %d %d %d\n' % (key[0], key[1], f, (A & B).length, (A | B).length, len_A, len_B))
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
