#!/usr/bin/env python
'''
============================================================
Main phasing program of a genotype data set for which
pedigree information is available.

It also supports imputation of the remaining missing
genotypes using genotype frequencies estimated from the data.

Created on July 26, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, time, util
from optparse import OptionParser
from impute.phasing.phase_trivial import trivial_phaser
from impute.phasing.phase_family import family_child_comparison_phaser, family_phaser
from impute.data import io
from impute.tools.param import PhaseParam
from impute.phasing.phase_core import new_phaser_chain
from impute.phasing.pre_processing import prepare_phaser
from impute.phasing.post_processing import save_stats_processor, impute_processor, fill_missing_processor
from impute.phasing import phase_distant

#---------------------------------------------
# Constants
#---------------------------------------------        
# Enumerated type for the fill-missing CLI option
class IMPUTE_OPTION: NONE, IMPUTE, IMPUTE_AND_FILL = range(3) 

# Number of phasing stages
NUM_STAGES = 6

#---------------------------------------------
# Methods
#---------------------------------------------           
def build_phasing_pipeline(options):
    '''Bulid and return the entire phasing processing chain by options - factory method.'''
    chain = __build_phaser(options) + __build_post_chain(options)
    return new_phaser_chain(chain, debug=options.debug, print_times=options.print_times)

def run_phasing_chain(phaser, problem, params=None):
    '''The main call that runs the phasing, stats saving , and post-processing as one long pipeline.
    Returns the populated request object.'''
    request = util.Struct(problem=problem, params=params if params else PhaseParam(),
                          g_orig=problem.genotype.data.copy(), stats=util.Struct())
    # Run phasing processing chain
    start = time.time()
    phaser.handle(request)
    t = time.time() - start
    request.stats.time = t
    
    return request 

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    prog = os.path.basename(sys.argv[0])
    usage = 'Usage: %s -i in_file -o out_file [-s stage] [-v]\n\n' \
        'Phase a genotype data set with pedigree information.\n\n' \
        'Type ''%s -h'' to display full help.' % (prog, prog)

    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--input', default=None,
                      help='Input file, NPZ format')
    parser.add_option('-o', '--output', default=None,
                      help='Output file, PLINK/NPZ format (NPZ assumed for .npz extension)')
    parser.add_option('-s', '--stage', default=0, type='int',
                      help='Run only this phasing stage (0=all; 1=trivial; 2=family: parent-offspring; \
                      3=family: children comparison)')
    parser.add_option('-f', '--impute', type='int', dest='impute', default=IMPUTE_OPTION.NONE,
                      help='Post-processing: do nothing (0), impute genotypes from called haplotypes \
                      (1), or impute and fill missing genotypes randomly from estimated frequencies (2)\
                      (default)')
    parser.add_option('-m', action='store_true', dest='min_output', default=False,
                      help=' Minimize output size: save only essential Problem fields in resulting npz')
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-p', '--processes', type='int', dest='num_processes', default=1,
                      help='Number of processes to spawn')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) > 0:
        print usage
        sys.exit(1)
    valid, options = __validate_options(options)
    if not valid:
        print usage
        sys.exit(1)
    return options
    
def __validate_options(options):
    '''Validate input options. Set default/overrides.'''
    # Argument validation
    if options.input is None or options.output is None:
        return False, options
    return True, options

def __load_problem(options):
    '''Load data from input files into a Problem object.'''
    # args = dict(prefix=options.prefix, tped=options.tped, tfam=options.tfam)
    # print args
    if options.debug: print 'Input file : %s ...\nOutput file: %s ...' % (options.input, options.output,)
    problem = io.read_npz(options.input)
    if options.debug: print 'Loaded', problem
    return problem
    
def __build_pre_phaser(options):
#    if options.stage != 0:
#        return []
    return [prepare_phaser()]

def __build_phaser(options):
    '''Bulid and return the phasing processing chain by options - factory method.'''
    if options.stage == 0:
        # 07-DEC-12: disabling stage 4 (family_sib_comparison_phaser) since it degrades phasing & imputation
        chain = [
                 prepare_phaser(),
                 trivial_phaser(),
                 family_phaser(),
                 family_child_comparison_phaser(),
                 phase_distant.family_sib_comparison_phaser(),
                 phase_distant.distant_phaser()
                 ]
    elif options.stage == 1: chain = [prepare_phaser(), trivial_phaser()]
    elif options.stage == 2: chain = [family_phaser()]
    elif options.stage == 3: chain = [family_child_comparison_phaser()]
    elif options.stage == 4: chain = [phase_distant.family_sib_comparison_phaser()]
    elif options.stage == 5: chain = [phase_distant.distant_phaser()]
    return chain

def __build_post_chain(options):
    '''Build optional post-processing.'''
#    if options.stage != 0:
#        return []
    chain = [save_stats_processor()]
    if options.impute >= IMPUTE_OPTION.IMPUTE: chain.append(impute_processor())
    if options.impute >= IMPUTE_OPTION.IMPUTE_AND_FILL: chain.append(fill_missing_processor())
    return chain

def __main(options):
    '''
    --------------------------------------------------
    Main program - accepts an options struct.
    --------------------------------------------------
    '''
    if options.debug: print 'Input options', options
    print 'Building phaser (stage = %d) ...' % (options.stage,)
    phaser = build_phasing_pipeline(options)
    
    if options.debug: print 'Reading data ...'
    problem = __load_problem(options)

    if options.debug: print 'Phasing ...'
    params = PhaseParam()
    params.update_from_struct(options)
    request = run_phasing_chain(phaser, problem, params)
    
    print ''
    request.stats.pprint()
    print ''

    if options.output is not None:
        if options.min_output:
            print 'Minimizing output size...'
            io.slim(problem)
        out_prefix, ext = os.path.splitext(options.output)
        if ext == '.npz':
            print 'Writing haplotype result to %s in NPZ format ...' % (options.output,)
            io.write_npz(problem, options.output)
            output_info = out_prefix + '.info.npz'
            print 'Writing problem info result to %s in NPZ format ...' % (output_info,)
            io.write_info_npz(problem.info, output_info)
        else:
            print 'Writing haplotype result to %s in PLINK format ...' % (options.output,)
            io.write_plink(problem, options.output, verbose=options.debug)
    return problem

#---------------------------------------------
# Main Program
#---------------------------------------------           
def main(**kwargs):
    '''Main program - accepts argument dictionary.'''
    # Default options
    options = util.Struct(pedigree=None, prefix=None, tped=None, tfam=None, input=None,
                          output=None, debug=False, stage=0, print_times=True,
                          impute=IMPUTE_OPTION.NONE, min_output=False, selected_samples=None)
    # Override with passed arguments
    options.update(**kwargs)
    valid, options = __validate_options(options)
    if not valid: raise ValueError('Bad options')
    return __main(options)
    
if __name__ == '__main__':
    '''Main program - accepts CLI arguments.'''
    __main(__parse_command_line_args())
