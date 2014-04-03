#!/usr/bin/env python
'''
============================================================
Convert our haplotype segment to IBD >= 1 segments between
the corresponding samples.

Created on February 12, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv, optparse, util, traceback, itertools, impute as im

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
    usage = 'Usage: %s [-v] <problem-file> [file name]\n\n' \
        'Convert our haplotype segment to IBD >= 1 segments between the corresponding\n' \
        'samples. Read from stdin, write to stdout.\n\n' \
        'Type ''%s -h'' to display full help.' % (prog, prog)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    options, args = parser.parse_args(sys.argv[1:])

    # Argument validation
    if len(args) < 1 or len(args) > 2:
        print usage
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
        
    # Read file name
    options.problem_file = args[0]
    options.file = sys.stdin if len(args) < 2 else open(args[1], 'rb')
    return options

read_lines = lambda in_file: csv.reader(in_file, delimiter=' ', skipinitialspace=True)
to_int_lines = lambda lines: ([int(y) for y in x] for x in lines)
groupby_hap_pair = lambda lines: ((k, list(tuple(x[2:4]) for x in g)) for k, g in itertools.groupby(lines, key=lambda x: x[4:8]))
sort_pair = lambda x: (min(x), max(x))
groupby_sample_pair = lambda hap_groups: ((k, list(map(lambda x: x[1], g))) 
                                          for k, g in itertools.groupby(hap_groups, key=lambda x: sort_pair((x[0][0], x[0][2]))))

union_segment_list = lambda lst: reduce(lambda A, B: im.segment.segment_set_op(A, B, 'or'), lst, [])
sample_pair_ibd = lambda sample_groups: ((k, im.segment.flatten(union_segment_list(g))) for k, g in sample_groups)

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    options = parse_command_line_args()
    try:
        problem = im.io.read_npz(options.problem_file)
        sample_id = problem.pedigree.sample_id
        a = read_lines(options.file)
        b = to_int_lines(a)
        c = groupby_hap_pair(b)
        d = groupby_sample_pair(c)
        e = sample_pair_ibd(d)
        for x in e:
            # IDs might make index pairs that were previously sorted unsorted => resort
            id1, id2 = sort_pair((sample_id[x[0][0]], sample_id[x[0][1]]))
            sys.stdout.write('%d %d %s\n' % (id1, id2, ' '.join('%d' % (y,) for y in x[1])))
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
