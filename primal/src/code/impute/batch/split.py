#!/usr/bin/env python
'''
============================================================
Split a PLINK BED data set into chromosome files.
Recode alleles so that 1=minor, 2=major.

Prerequisites: PLINK must be installed and be on the PATH.

Created on January 13, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, traceback, numpy as np, util, impute as im, impute.batch.batch_util as bu, db_gene
from optparse import OptionParser

#---------------------------------------------
# Constants
#---------------------------------------------
# Program name
PROGRAM = os.path.basename(sys.argv[0])

#---------------------------------------------
# Methods
#---------------------------------------------

#---------------------------------------------
# Private Methods
#---------------------------------------------

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    usage = 'Usage: %s <plink_set> <chr> <out>\n\n' \
        'Read plink_set.bim file and compute part size per chromosome. Save in plink_set.part;\n' \
        'or split a PLINK BED file into chromosome parts accordingly.\n' \
        'part_size is the slice size in Mb.\n\n' \
        'Example 1: Print part counts by chromosome: %s genome 22 genome-chr22\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-v', '--debug'        , action='store_true'  , dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-r', '--recode'      , action='store_true'  , dest='recode', default=True,
                      help='Recode alleles to 1=minor, 2=major (if False, a random assignment to 1,2 is made)')
    parser.add_option('-i', '--id-coef', type=str, dest='id_coef', default=None,
                      help='Identity coefficient file for all sample pairs. Format: id1 id2 lam delta1...delta9. If empty, defaults to plink_set.id')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 3:
        print usage
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
    input, chrom, out = args[0], int(args[1]), args[2]  # @ReservedAssignment
    if chrom < 1 or chrom > 22:
        print usage
        print('\nMust specify a chromosome number in 1..22.')
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
        
    try:
        util.mkdir_if_not_exists(os.path.dirname(out))
        # Use [cM] as genetic distance unit
        plink_cmd_base = '%s --bfile %s --chr %d --out %s' % (bu.PLINK, input, chrom, out)
        
        if options.recode:
            # First, compute allele frequencies with PLINK  
            util.run_command(plink_cmd_base + ' --nonfounders --freq')
            # Convert frequencies file that to a reference allele recoding
            # file (a file containing the list of SNPs and their minor allele letter)
            bu.frq_to_minor_file(out + '.frq', out + '.mnr') 
            # Finally, convert binary PLINK to a 12-recoded TPED, where 1=minor allele for each SNP                 
            util.run_command('%s --transpose --recode12 --reference-allele %s.mnr' % (plink_cmd_base, out))
        else:
            # No recoding, just convert binary to 2-recoded TPED. PLINK assigns "1" to
            # the first-encountered allele in the file for each SNP.
            util.run_command('%s --transpose --recode12' % (plink_cmd_base,))
            
        # Extract chromosome frames, save to frm file
        frames = db_gene.snp.ld_graph.read_frames(input + '.frm')
        db_gene.snp.ld_graph.write_frames(db_gene.snp.ld_graph.Frames((chrom, x) for x in frames[chrom]), out + '.frm')

        # Calculate lambda(f) and save to lam file
        lam = im.hap_lambda.lambda_vs_f(options.id_coef if options.id_coef else (input + '.id'))
        F, L, _ = im.hap_lambda.lambda_mean(lam)
        np.savetxt(out + '.lam', np.array([F, L]).transpose())
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
