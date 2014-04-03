#!/usr/bin/env python
'''
============================================================
Convert a Problem object from PLINK to NPZ format or vice
versa.

Created on August 5, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, util
from optparse import OptionParser
from impute.data import io

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    prog = os.path.basename(sys.argv[0])
    usage = 'Usage: %s -t <npz|plink> [-p <pedigree_file>] -g <genotype_file_prefix> -z <npz_file> [-v]\n\n' \
        'Convert a PLINK file set with (pedigree_file, file_prefix{.tped,.tfam,.hap.tped,.info})\n' \
        'to an npz_file (if -t npz is specified) or vice versa (if -t plink is specified). If\n' \
        'no .hap.tped file is found, an empty haplotype set is assumed.\n' \
        '-p is needed only with ''-t npz''.\n' \
        'no .info file is found, an empty ProblemInfo object is assumed.\n\n' \
        'Type ''%s -h'' to display full help.' % (prog, prog)

    parser = OptionParser(usage=usage)
    parser.add_option('-p', '--pedigree', default=None,
                      help='Pedigree file (typed+untyped)')
    parser.add_option('-q', '--pedigree-genotyped',
                      type='str', dest='pedigree_genotyped', default=None,
                      help='Pedigree of genotyped samples')
    parser.add_option('-g', '--prefix', default=None,
                      help='Genotype data set prefix (plink .tfam,.tped required; .hap.tped optional)')
    parser.add_option('-z', '--npz', default=None,
                      help='NPZ file name')
    parser.add_option('-t', '--target', default=None,
                      help='Target format (npz/plink)')
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) > 0:
        print usage
        sys.exit(1)
    valid, options, message = __validate_options(options)
    if not valid:
        print usage
        print message
        sys.exit(1)
    return options
    
def __validate_options(options):
    '''Validate input options. Set default/overrides.'''
    # Set defaults/overrides
    if options.prefix is None or options.npz is None or options.target is None:
        return (False, options, 'Must specify prefix, npz file or target')
    if options.target == 'npz' and options.pedigree is None:
        return (False, options, 'Must have a pedigree if converting from plink')
    options.tfam = options.prefix + '.tfam'
    options.tped = options.prefix + '.tped'
    options.haplotype = options.prefix + '.hap.tped'
    options.info = options.prefix + '.info'
    options.frames = options.prefix + '.frm'
    return (True, options, None)
   
def __main(options):
    '''
    --------------------------------------------------
    Main program - accepts an options struct.
    --------------------------------------------------
    '''
    if options.target == 'npz':
        haplotype = options.prefix + '.hap.tped'
        haplotype = haplotype if os.path.isfile(haplotype) else None 
        return io.plink_to_npz(options.prefix, options.npz, options.pedigree, haplotype=haplotype,
                               info=None,verbose=options.debug,
                               pedigree_genotyped=options.pedigree_genotyped)
    else:
        return io.npz_to_plink(options.npz, options.prefix, verbose=options.debug)

#---------------------------------------------
# Main Program
#---------------------------------------------
def main(**kwargs):
    '''Main program - accepts argument dictionary.'''
    # Default options
    options = util.Struct(pedigree=None, prefix=None, tped=None, tfam=None, haplotype=None,
                          out=None, debug=False, target=None, pedigree_genotyped=None)
    # Override with passed arguments
    options.update(**kwargs)
    valid, options, error_msg = __validate_options(options)
    if not valid: raise ValueError('Bad options: %s' % (error_msg,))
    return __main(options)
    
if __name__ == '__main__':
    '''Main program - accepts CLI arguments.'''
    __main(__parse_command_line_args())
