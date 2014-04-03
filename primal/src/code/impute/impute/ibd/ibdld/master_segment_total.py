#!/usr/bin/env python
'''
============================================================
Calculate the total % of IBD between each individual pair
in a master segment file.

Created on February 13, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv, optparse, impute as im, util, traceback
from impute.ibd.distant.ibd_hmm_hap import kinship, p_ibd

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
    usage = 'Usage: %s [-v]\n\n' \
    'Calculate the total %% of IBD between each individual pair in a master segment file.\n' \
    'Reads from standard input; writes to standard output.\n\n' \
    'Type ''%s -h'' to display full help.' % (prog, prog)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-i', '--idcoef-file', type='str'  , dest='id_coef_file',
                      default=None,
                      help='Path to identity coefficient file (if None, using default location)')
    options, args = parser.parse_args(sys.argv[1:])
    
    # Argument validation
    if len(args) != 0:
        print usage
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
        
    # Read file name
    #options.file = args[0]
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
    param = im.PhaseParam()
    if options.id_coef_file:
        param.id_coef_file = options.id_coef_file
    try:
        for line in csv.reader(sys.stdin, delimiter=' ', skipinitialspace=True):
            key, A = parse_line(line)
            # Output statistics. Only log pairs for which samples were found
            len_A = A.length
            if len_A != 0:
                #f = param.kinship(key[0], key[1])
                _, Delta = param.id_coefs(key[0], key[1])
                f = kinship(Delta)
                pibd = p_ibd(Delta)
                sys.stdout.write('%d %d %f %f %d\n' % (key[0], key[1], f, pibd, len_A))
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
