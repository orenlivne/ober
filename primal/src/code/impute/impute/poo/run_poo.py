#!/usr/bin/env python
'''
============================================================
Main driver for calculating parent-of-origin phases.

Created on August 8, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, impute as im, numpy as np, optparse, util

#---------------------------------------------
# Methods
#---------------------------------------------           
def poo_alignment(problem, params):
    '''Phase all homozygous SNPs: find all homozygous SNPs in the genotype set
    and set their corresponding haplotypes to the genotypes.''' 
    h = problem.haplotype
    ped = problem.pedigree
    if params.debug:
        print 'Running POO alignment, sample threshold %.2f, snp threshold %.2f, SNP step size %d' % \
        (params.poo_threshold, params.poo_snp_measure_threshold, params.poo_snp_step_size)
    # Calculate POO phases 
    m = im.poo.determine_poo(problem.info.snp['chrom'][0], params=params)
    # Save POO phases of all samples.
    # Non-quasi-founder samples: phase = 1 (no hap flipping required)
    # Quasi-founder (QF) phase: 1 if not need to flip, -1 if need to flip, 0 if undetermined
    poo_phase = np.zeros((problem.num_samples,), dtype=np.byte)
    non_founders = np.where([any((y < ped.num_genotyped) for y in ped.graph.predecessors(x)) for x in xrange(ped.num_genotyped)])[0]
    ind = (m[:, 1] >= params.poo_threshold)
    aligned_qf = m[ind, 0].astype(np.uint)
    poo_phase[non_founders] = 1
    poo_phase[aligned_qf] = m[ind, 6].astype(np.byte)
    h.poo_phase = poo_phase
    if params.debug:
        print 'aligned_quasi_founders %d non_founders %d total aligned %d of %d' % \
        (len(np.where(ind)[0]), len(non_founders), len(h.aligned_samples), ped.num_genotyped)
    # Set parent-of-origin genotype tags
    h.hap_type[:] = im.constants.PHASED
    h.hap_type[:, h.aligned_samples] = im.constants.PHASED_WITH_ORIGIN

    return False

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __determine_poo_sample((pedigree, params, sample, lock)):
    # m = a.flip_measure_sample(sample, threshold=params.poo_snp_measure_threshold, snp_step_size=params.poo_snp_step_size, debug=params.debug)[0]
    m = (params, sample)
    util.write_with_lock('%s\n' % (repr(m),), lock)
    return m

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __parse_command_line_args(argv):
    '''Parse and validate command-line arguments.'''
    PROGRAM = os.path.basename(argv[0])
    usage = 'Usage: %s [flags] <in.npz> <out.npz>\n\n' \
        'Load a problem from in.npz, align POO phases and save them with the problem in out.npz.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-d', '--debug', type='int', dest='debug', default=0,
                      help='Debug Level (0=quiet; 1=summary; 2=full debug)')
    parser.add_option('-i', '--ibd-index', type='str'  , dest='index_file',
                      default=os.environ['OBER_OUT'] + '/index_segments',
                      help='Path to IBD index (top directory)')
    parser.add_option('-p', '--processes', type='int', dest='num_processes', default=1,
                      help='Number of processes to spawn')
    parser.add_option('-s', '--snp-step-size', type='int', dest='poo_snp_step_size', default=1,
                      help='POO SNP downsampling step size')
    options, args = parser.parse_args(argv[1:])
    if len(args) != 2:
        print usage
        sys.exit(1)
    options.input = args[0] 
    options.output = args[1]
    return args, options

#---------------------------------------------
# Main Program
#---------------------------------------------           
if __name__ == '__main__':
    args, options = __parse_command_line_args(sys.argv)
    params = im.PhaseParam()
    params.update_from_struct(options)
    
    problem = im.io.read_npz(options.input)
    poo_alignment(problem, params)
    im.io.write_npz(problem, options.output)
