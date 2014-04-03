#!/usr/bin/env python
'''
============================================================
Read a PLINK BIM file and the SNP database, and export text
files containing all information required for phasing that
is not in the PLINK data set. 
 
Created on January 21, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, traceback, numpy as np, util, impute.batch.batch_util as bu, db_gene
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
    usage = 'Usage: %s <bim-file>\n\n' \
        'Read a BIM file and extract information into phasing input files.\n' \
        '- The BIM file is copied to a new one with a populated genetic coordinate column;\n' \
        '- An .frm file is also created.\n' \
        '- The new BIM is called bim-file.new.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-v', '--debug'        , action='store_true'  , dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-d', '--db-url'       , type='str'           , dest='db_url', default=db_gene.DEFAULT_URL,
                      help='SNP database URL')
    parser.add_option('-o', '--out'          , type='str'           , dest='out_base_name', default=None,
                      help='Output PLINK data set base name')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 1:
        print usage
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
    input_file = args[0]
    if not options.out_base_name: options.out_base_name = os.path.splitext(input_file)[0]
            
    try:
        # Initialize
        daos = db_gene.snp.snp_db_dao.Daos(url=options.db_url)
        util.mkdir_if_not_exists(os.path.dirname(options.out_base_name))
        
        # Set genetic distance column in BIM file (read locations from snp db) and save a new copy of it
        snp_data = np.genfromtxt(input_file,
                                 dtype=[
                                        ('chrom', np.uint8),  # Chromosome # containing the SNP
                                        ('name', np.chararray),  # SNP name (e.g., 'rs...')
                                        ('dist_cm', np.float),  # Genetic position [CENTI-Morgans!!]
                                        ('base_pair', np.uint),  # Base pair position on chromosome
                                        ('allele1', np.chararray),
                                        ('allele2', np.chararray)
                                        ])
        snp_names = snp_data['name']
        a = dict((x.name, x) for x in daos.snp_dao.get_snps_iter(snp_names))
        # Note: our genetic distance unit is cM 
        snp_data['dist_cm'] = map(lambda x: x if x else 0.0, ((a[x].genetic_pos if a.has_key(x) else None) for x in snp_names))
        np.savetxt(options.out_base_name + '.bim.new', snp_data, fmt='%d\t%s\t%f\t%d\t%s\t%s')

        # For each chromosome in the PLINK file: load LD data, generate frame numbers, save them to a file
        (c, s), frames = bu.get_bim_metadata(open(input_file, 'rb')), util.mdict()
        for chrom in s.iterkeys():
            for x in db_gene.snp.ld_graph.frames(chrom, s[chrom], daos.snp_dao, daos.ld_dao):
                frames[chrom] = x
        with open(options.out_base_name + '.frm', 'wb') as frm_file:
            db_gene.snp.ld_graph.write_frames(frames, frm_file)
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
