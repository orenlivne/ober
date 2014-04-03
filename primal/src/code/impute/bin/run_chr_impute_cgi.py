#!/usr/bin/env python
'''
============================================================
Impute Hutterites CGI WGS data at SNP locations on a
single Chromosome. Read from file; write to standard output
in CGI format.

Created on February 15, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv, util, db_gene, numpy as np, traceback, impute as im, optparse, itemutil

#---------------------------------------------
# Constants
#---------------------------------------------

#---------------------------------------------
# Methods
#---------------------------------------------
def __parse_command_line_args(argv):
    '''Parse and validate command-line arguments.'''
    PROGRAM = os.path.basename(argv[0])
    usage = 'Usage: %s [flags] <genotype-file> <id-file> <phasing-problem-file.npz> <ibd-segment-file>\n\n' \
        'Impute Hutterites CGI whole-genome sequencing data at SNP locations on a\n' \
        'single chromosome. Read genotype and CGI sample ID files.\n'\
        'Write to standard output in CGI format.\n' \
        'If genotype-file = ''-'', takes input from stdin.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-a', '--algorithm', type='str', dest='algorithm', default='index',
                      help='Imputation algorithm: interval_tree|index.')
    parser.add_option('-b', '--block-size', type='int', dest='block_size', default=3,
                      help='SNP block size to be bulk-read from the data and imputed in one call')
    parser.add_option('-c', '--chr'          , type='int'           , dest='chrom', default=0,
                      help='Chromosome number to process. Overrides other options if non-zero. 0=do chromosomes in the range specified by the start, end chromosomes')
    parser.add_option('-d', '--debug', type='int', dest='debug', default=0,
                      help='Debug Level (0=quiet; 1=summary; 2=full debug)')
    parser.add_option('-i', '--index-file', type='str'  , dest='index_file',
                      default=os.environ['OBER_DATA'] + '/cgi/README.assembly_sample_subject.csv',
                      help='CGI+FINDIV ID index file')
    parser.add_option('-o', '--output-id', action='store_true'  , dest='out_id',
                      default=False,
                      help='If set, outputs training set IDs instead of imputing. In that case, only one mandatory argument is expected: phasing-problem-file.npz')
    parser.add_option('-s', '--debug-sample', type='int', dest='debug_sample', default= -1,
                      help='Debug at sample %d (-1=no sample)')
    parser.add_option('-j', '--inject-ids', action='store_true'  , dest='inject_ids', default=False,
                      help='If true, the IDs read from id-file are interpreted as sample FINDIVs rather than' \
                      'CGI GSM sequence identifiers')
    options, args = parser.parse_args(argv[1:]) 
    if len(args) != (1 if options.out_id else 4):
        print usage
        sys.exit(1)
    if options.chrom < 1 or options.chrom > 22:
        print usage
        print('\nMust specify a chromosome number in 1..22.')
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
    return args, options

'''Read CGI sample IDs and convert to our sample IDs.'''
__read_sample_ids = lambda cgi_id_file, cgi_id_to_sample_id: util.sort_with_index(np.array([cgi_id_to_sample_id[cgi_id] for cgi_id in __csv_reader(cgi_id_file).next()[8:]]))
__read_sample_ids_inject_ids = lambda cgi_id_file, cgi_id_to_sample_id: util.sort_with_index(np.loadtxt(cgi_id_file, dtype=int))

'''Parse a CSV file name f.'''
__csv_reader = lambda f: csv.reader(open(f, 'rb') if isinstance(f, str) else f, skipinitialspace=True, delimiter='\t')

def __imputer(problem, segment_location, options):
    '''A factory method that selects that imputation algorithm implementation.'''
    algorithm = options.algorithm
    if algorithm == 'interval_tree': return _ImputerIntervalTree(problem, segment_location, options)
    elif algorithm == 'index': return _ImputerIndex(problem, segment_location, options)
    else: raise ValueError('Unsupported imputation algorithm ''%s''' % (algorithm,))

####################################################################################
class _Imputer(object):
    '''Base class for imputation implementations.'''
    
    def __init__(self, problem, options):
        '''Read phased haplotypes and the IBD segment dictionary.'''
        self._haplotype = problem.haplotype  # May or may not be used by sub-class implementation
        self._debug = options.debug
        self._debug_sample = options.debug_sample 
    
    def run(self, t):
        '''Impute a f SNPs using the phased haplotypes in self._haplotype and CGI
        genotype data in t. t.imputed_data is updated in-place with the imputation results.'''
        raise ValueError('Must be implemented by sub-classes')

####################################################################################
class _ImputerIntervalTree(_Imputer):
    '''Version 1 (slow): imputes using interval tree queries of a SmartSegmentSet.'''

    def __init__(self, problem, segment_file, options):
        '''Read phased haplotypes and the IBD segment dictionary.'''
        super(_ImputerIntervalTree, self).__init__(problem, options)
        self._ibd = im.smart_segment_set.SmartSegmentSet.load(problem.pedigree.num_genotyped, segment_file)
    
    def run(self, t):
        '''Impute a f SNPs using the phased haplotypes in self._haplotype and CGI
        genotype data in t. t.imputed_data is updated in-place with the imputation results.'''
        im.imputation.iibd.impute(self._haplotype, self._ibd, t, debug=self._debug)

####################################################################################
class _ImputerIndex(_Imputer):
    '''Version 1 (slow): imputes using interval tree queries of a SmartSegmentSet.'''
    def __init__(self, problem, segment_location, options):
        '''Read phased haplotypes and the IBD segment dictionary.'''
        super(_ImputerIndex, self).__init__(problem, options)
        self._ibd = im.index.segment_index.SegmentIndex(segment_location, chrom=options.chrom)
    
    def run(self, t):
        '''Impute a f SNPs using the phased haplotypes in self._haplotype and CGI
        genotype data in t. t.imputed_data is updated in-place with the imputation results.'''
        im.imputation.impute_ibd_index.impute(self._ibd, t, debug=self._debug, debug_sample=self._debug_sample)

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    args, options = __parse_command_line_args(sys.argv)
    if options.debug:
        print options
        print args   
    if options.out_id:
        problem_file = args[0]
    else:
        input_file = sys.stdin if args[0] == '-' else open(args[0], 'rb')
        cgi_id_file, problem_file, segment_location = args[1:4]
    try:
        # Convert sample IDs to sample indices 
        problem = im.io.read_npz(problem_file)
        if options.out_id:
            #--------------
            # SAVE INDEX
            #--------------
            np.savetxt(sys.stdout, problem.pedigree.sample_id, '%d')
        else:
            #--------------
            # IMPUTE
            #--------------
            if options.debug:
                sys.stdout.write('Imputing chromosome %d\n' % (options.chrom,))
            # Create a CGI-ID-to-our-sample-ID (FINDIV) dictionary
            cgi_id_to_sample_id = db_gene.cgi.ids.cgi_id_to_sample_id(options.index_file)
            sample_id_reader = __read_sample_ids_inject_ids if options.inject_ids else __read_sample_ids
            sample_id, index = sample_id_reader(cgi_id_file, cgi_id_to_sample_id)
            training_id = cgi_id_to_sample_id.values()
            training_index = np.array([problem.pedigree.node_of[x] for x in training_id])
            if options.debug >= 1:
                print 'sample_id', sample_id
                print 'index', index
                print 'training_id', training_id
                print 'training_index', training_index
                j = np.where(np.array(training_id) == 162152)[0]
                print j
                print training_id[j], training_index[j]
    
            imputer = __imputer(problem, segment_location, options)
    
            # Loop over blocks of SNPs and impute each block        
            # TODO: optimize block size - maybe replace by a smart iterator that reads blocks of 
            # SNPs that lie between two Affy SNPs to save IBD segment query time (one query for all of them)
            # Less of a problem in imputation V2. Can take larger block sizes.  
            for block in itemutil.group_iterates(__csv_reader(input_file), options.block_size):
                # Read genotype data into Genotype object
                g = im.cgi.io_cgi.read_genotype(options.chrom, sample_id, index, block)
                t = im.imputation.ImputationSet(problem.pedigree, g)
                if options.debug >= 1:
                    print g
                    print t
                imputer.run(t)
                # Write imputed genotypes  in CGI format 
                im.cgi.io_cgi.write_imputed(t, sys.stdout, debug=options.debug >= 3, poo_phase=problem.haplotype.poo_phase)
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
