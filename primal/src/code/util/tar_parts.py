#!/usr/bin/env python
'''
============================================================
Archive all part files into a single output file.

Created on Septmeber 21, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, zipfile, datetime
from optparse import OptionParser

#---------------------------------------------
# Constants
#---------------------------------------------
# Program name
PROGRAM = os.path.basename(sys.argv[0])

def print_zip_info(archive_name):
    zf = zipfile.ZipFile(archive_name)
    for info in zf.infolist():
        print info.filename
        print '\tComment:\t', info.comment
        print '\tModified:\t', datetime.datetime(*info.date_time)
        print '\tSystem:\t\t', info.create_system, '(0 = Windows, 3 = Unix)'
        print '\tZIP version:\t', info.create_version
        print '\tCompressed:\t', info.compress_size, 'bytes'
        print '\tUncompressed:\t', info.file_size, 'bytes'
        print
    
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse command-line arguments
    usage = 'Usage: %s <infile_basename> <part_count_file> <chr> <out_file>\n\n' \
        'Reduces a file set of averages into a global average.\n' \
        'Type ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-v', '--debug'        , action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    (options, args) = parser.parse_args(sys.argv[1:])
    if len(args) != 4:
        print usage
        sys.exit(1)
    (in_file, part_count_file, chrom, out_file) = args
    chrom = int(chrom)

    print_zip_info(in_file, part_count_file, chrom, out_file)
    if options.debug:
        print_zip_info(out_file)
