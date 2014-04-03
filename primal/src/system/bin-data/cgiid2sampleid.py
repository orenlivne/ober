#!/usr/bin/env python
'''
============================================================
Convert Hutterites WGS CGI data at SNP locations read from
standard input to PLINK TPED format.

Prerequisites: tabix

Created on November 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv, db_gene
from optparse import OptionParser

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <sample ID index file>\n' \
        'Convert CGI ID to sample ID.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 1:
        print usage
        sys.exit(1)
    
    # Read CGI-ID-to-our-sample-ID (FINDIV) dictionary
    sample_index_file = os.environ['OBER'] + '/data/cgi/README.assembly_sample_subject.csv' \
    if len(args) < 1 else args[0]
    cgi_id_to_sample_id = db_gene.cgi.ids.cgi_id_to_sample_id(sample_index_file)
    
    for line in csv.reader(sys.stdin, skipinitialspace=True):
        print line
        if line:
            print cgi_id_to_sample_id[line]
