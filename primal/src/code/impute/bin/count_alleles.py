#!/usr/bin/env python
'''
============================================================
Calculate allele frequencies in a .gxn (Gaixin-format) file. 

GXN format: snp name, sample ID, allele1, allele2.
Typically, the 1 allele is the minor and 2 is the major.

Row order does NOT matter.

SNP_1 SAMPLE_ID_1 A1 A2
...
SNP_1 SAMPLE_ID_N A1 A2
SNP_2 SAMPLE_ID_1 A1 A2
...
SNP_2 SAMPLE_ID_N A1 A2
...
SNP_M SAMPLE_ID_1 A1 A2
...
SNP_M SAMPLE_ID_N A1 A2

Created on October 18, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv
from optparse import OptionParser

#---------------------------------------------
# Constants
#---------------------------------------------
# Allele indices in count dictionary entries
A1, A2 = range(2)

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s\n' \
        'Calculate allele frequencies of each SNP in a gxn-format standard input.\n' \
        'For each SNP, lists: SNP name, A1 count, A2 count, A1 frequency, A2 frequency.\n' \
        'List all SNPs: cat hutt_phased_chr22_part0.gxn | count_alleles.py\n' \
        'List only SNPs where the 1-allele is not minor: cat hutt_phased_chr22_part0.gxn | count_alleles.py -m\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-m', '--flag-maf-not-one', action='store_true'  , dest='flag_maf_not_one', default=False,
                      help='If true, list only the SNPs in which the 1-allele is not the minor one')
    (options, args) = parser.parse_args(sys.argv[1:])
    if len(args) != 0:
        print usage
        sys.exit(1)
    
    count = {}
    for line in csv.reader(sys.stdin, delimiter=' ', skipinitialspace=True):
        if line:
            (snp, a1, a2) = (line[0], int(line[1]), int(line[2]))
            count.setdefault(snp, [0, 0])
            count[snp][A1 if a1 == 1 else A2] += 1
            count[snp][A1 if a2 == 1 else A2] += 1

    try:
        for (k, v) in sorted(count.items(), key=lambda x: x[0]):
            total = 1.0 * (v[0] + v[1])
            (f1, f2) = (v[0] / total, v[1] / total)
            if f1 >= 0.5 or not options.flag_maf_not_one: 
                sys.stdout.write('%-20s %4d %4d %7.5f %7.5f\n' % (k, v[0], v[1], v[0] / total, v[1] / total))
        sys.stdout.flush()
    except (IOError, OSError):
        sys.exit(141)