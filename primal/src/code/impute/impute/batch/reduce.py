#!/usr/bin/env python
'''
============================================================
Reduces all chromosome part data sets into a single data set,
or all chromosomes into the final data set.

Created on Septmeber 25, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, traceback, shutil, tempfile, numpy as np, impute.batch.batch_util as bu, util
from optparse import OptionParser
from util import mkdir_if_not_exists

#---------------------------------------------
# Constants
#---------------------------------------------
# Program name
PROGRAM = os.path.basename(sys.argv[0])
# PLINK binary data set file extensions
EXTENSIONS = ['bed', 'bim', 'fam']

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __parse_part_type(part_type_str):
    '''Parse the part type CLI argument into the corresponding enumerated type.'''
    if part_type_str == 'chr':
        return bu.CHROMOSOME
    elif part_type_str == 'part':
        return bu.PART
    else:
        raise ValueError('Unsupported partition type ''%s''' % (part_type_str,))

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse command-line arguments
    usage = 'Usage: %s -s <start_part> -e <stop_part> <plink-set> <part-type> <out-plink-set>\n\n' \
        'Reduce chromosome part PLINK binary sets to a single binary set.\n' \
        'Type ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-s', '--start-part'          , type='int', dest='start', default=None,
                      help='Starting part number (inclusive)')
    parser.add_option('-e', '--stop-part'          , type='int', dest='stop', default=None,
                      help='Ending part number (not inclusive)')
    (options, args) = parser.parse_args(sys.argv[1:])
    if len(args) != 3:
        print usage
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
    if options.start is None or options.stop is None:
        print 'Must specify start and stop'
        print usage
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
    (in_file, part_type, out_file) = args
    part_type = __parse_part_type(part_type)
    num_parts = options.stop - options.start
    mkdir_if_not_exists(os.path.dirname(out_file))
    
    try:
        # Merge PLINK data sets. If there's one part, nothing to merge, just copy the files.
        part_names = bu.partition_names(in_file, part_type,
                                        parts=xrange(options.start, options.stop)).values()
        first_part_name = part_names[0]
        print 'Reducing, num_parts', num_parts
        if num_parts == 1:
            for ext in EXTENSIONS:
                shutil.copy(first_part_name + '.' + ext, out_file + '.' + ext)
        else:
            # Prepare PLINK merge command input file
            f = tempfile.NamedTemporaryFile(delete=False)
            for name in part_names:
                for ext in EXTENSIONS:
                    f.write('%s.%s ' % (name, ext))
                f.write('\n')
            f.close()
            # Run PLINK merge of the first part with the rest
            cmd = '%s --bfile %s --merge-list %s --make-bed --out %s' % \
            (bu.PLINK, first_part_name, f.name, out_file)
            util.run_command(cmd)
            os.remove(f.name)
            
        # Merge stat npz files into a single npz file with an array of stats
        # TODO: convert this code to list comprehension : a tuple/dictionary of all fields over all files 
        (out_stats, out_info, out_pedigree) = ([], [], [])
        first = True
        for name in part_names:
            files = np.load(name + '.stats.npz')
            out_stats.append(files['stats'][0])
            out_info.append(files['info'][0])
            # Save only one copy of the pedigree (from the first part) in the merged file
            if first:
                out_pedigree = files['pedigree'][0]
                first = False
        np.savez(out_file + '.stats',
                 stats=np.array([out_stats]),
                 info=np.array([out_info]),
                 pedigree=np.array([out_pedigree]))
        
        # Create a file that tells observers that we're done
        if part_type == bu.CHROMOSOME:
            open(out_file + 'done', 'wb').close() 
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
