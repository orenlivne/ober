#!/usr/bin/env python
'''
============================================================
Phase a PLINK data set (typically, a chromosome part) and
imputed the genotypes from fully-called haplotypes.

tped -> -> npz -> phase -> npz -> {tped, stats npz} ->
{bed, stats npz}.


Created on September 24, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, traceback, numpy as np, impute.batch.batch_util as bu, util
from optparse import OptionParser
from impute.phasing import phase
from impute.data import io, io_genotype
from util import mkdir_if_not_exists

#---------------------------------------------
# Constants
#---------------------------------------------
# Program name
PROGRAM = os.path.basename(sys.argv[0])

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
    usage = 'Usage: %s <plink-set> <full-pedigree-tfam> <out-plink-set>\n\n' \
        'Phase a PLINK TPED data set. The full pedigree contains both genotyped\n' \
        'and non-genotyped individiduals.\n\n' \
        'Outputs a binary PLINK imputed genotype file and a phasing statistics NPZ file.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-v', '--debug'          , action='store_true'  , dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-t', '--stage'          , type='int'           , dest='stage', default=0,
                      help='Run only a this phasing stage')
    parser.add_option('-f', '--impute', type='int', dest='impute', default=phase.IMPUTE_OPTION.NONE,
                      help='Post-processing: do nothing (0), impute genotypes from called haplotypes \
                      (1), or impute and fill missing genotypes randomly from estimated frequencies (2)\
                      (default)')
    parser.add_option('-r', '--recode'      , action='store_true'  , dest='recode', default=False,
                      help='Recode alleles to 1=minor, 2=major (if False, allele coding is kept intact)')
    parser.add_option('-g', '--out-gxn'     , type='str'           , dest='out_gxn', default=bu.ARG_NONE,
                      help='Output directory of GXN files (if not specified, writes to same directory as out-plink-set''s)')
    (options, args) = parser.parse_args(sys.argv[1:])
    options.print_times = True
    if options.out_gxn.startswith(bu.ARG_NONE):
        options.out_gxn = None
    if len(args) != 3:
        print usage
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
    
    try:
        # Prepare file names, create directories
        (base_name, pedigree_file, out_base_name) = args
        mkdir_if_not_exists(os.path.dirname(out_base_name))
        if options.out_gxn:
            mkdir_if_not_exists(os.path.dirname(options.out_gxn))
        else:
            options.out_gxn = out_base_name
    
        npz_file = base_name + '.npz'
        
        # Convert plink tped -> npz
        problem = io.read_plink(prefix=base_name, pedigree=pedigree_file, haplotype=None,
                                verbose=options.debug)
    
        # Phase, impute, fill missing
        phaser = phase.build_phasing_pipeline(options)      
        request = phase.run_phasing_chain(phaser, problem)
        stats = request.stats
        print ''
        stats.pprint()
        print ''
    
        # Convert phased npz -> plink tped. Save only genotypes (haplotypes may need to be saved in the stats
        # object as a hash table for 'coloring the pedigree' at a later stage.
        genotype_file = out_base_name + '.tped'
        io.write_plink(problem, out_base_name, verbose=True,
                       save_node_type=False, save_genotype=True, save_haplotype=False, save_error=False)
        
        # Save statistics and phasing metadata in a separate npz
        np.savez(out_base_name + '.stats', 
                 stats=np.array([stats]), 
                 info=np.array([problem.info]), 
                 pedigree=np.array([problem.pedigree]))

        plink_cmd_base = '%s --tfile %s' % (bu.PLINK, out_base_name,)
        if options.recode:
            # First, compute allele frequencies with PLINK  
            util.run_command('%s --nonfounders --freq --out %s' % (plink_cmd_base, out_base_name))
            # Convert frequencies file that to a reference allele recoding
            # file (a file containing the list of SNPs and their minor allele letter)
            bu.frq_to_minor_file(out_base_name + '.frq', out_base_name + '.mnr') 

            # Then convert binary PLINK to a recoded 12-recoded TPED, where 1=minor allele for each SNP
            out_recoded = out_base_name + '.recoded'                 
            util.run_command('%s --transpose --recode12 --reference-allele %s.mnr --out %s' % \
                           (plink_cmd_base, out_base_name, out_recoded))

            # Reload the recoded problem
            for ext in ('nof', 'tped', 'tfam'):
                os.rename(out_recoded + '.' + ext, out_base_name + '.' + ext)
            genotype = io_genotype.read('plink', 'genotype', tped=out_base_name + '.tped', load_ids=False)
        else:
            genotype = problem.genotype
            
        # Write problem to file in our (npz)
        io.write_npz(problem, out_base_name + '.npz')
        # Write genotypes Gaixin formats; she uses those separate files
        io_genotype.write('gaixin', genotype, options.out_gxn + '.gxn',
                          sample_id=problem.pedigree.sample_id_genotyped)
             
        # Convert plink tped to bed; delete the tped set
        util.run_command('%s --make-bed --out %s' % (plink_cmd_base, out_base_name))
        for ext in ('nof', 'pdg.tfam', 'tped', 'tfam', 'info'):
            os.remove(out_base_name + '.' + ext)
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
