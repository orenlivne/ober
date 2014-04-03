#!/usr/bin/env python
'''
============================================================
Run phasing in stages on a single chromosome part BED file. 

Created on August 16, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, impute as im, itertools, csv, optparse, traceback, util, time, operator
from multiprocessing import Pool, Manager
from impute.ibd.distant.ibd_hmm_hap import IbdProblem
import itemutil

#---------------------------------------------
# Methods
#---------------------------------------------
'''Read haplotype pairs from file.'''
read_pairs = lambda input_file: (((i, a), (j, b)) for ((i, j), (a, b)) in 
                                 itertools.product(((int(line[0]), int(line[1])) for line in 
                                                    csv.reader(input_file, delimiter=' ', skipinitialspace=True) if line),
                                                   itertools.product(im.constants.ALLELES, im.constants.ALLELES))
                                 if (i, a) != (j, b))

#---------------------------------------------
# Private Methods
#---------------------------------------------
def process_pair((ibd_problems, options, params, lock)):
    '''Process sample pairs - outputs their IBD segments.'''
    processor = PairProcessor(options, params, lock)
    return processor.run(ibd_problems)

####################################################################################
class PairProcessor(object):
    '''Processes sample pairs - outputs their IBD segments.
    Load phased data - once per process. Not ideal but better than passing in this object,
    which would have been copied, which is also slow.'''

    #---------------------------------------------
    # Constants
    #---------------------------------------------
    '''Size of segment buffer.'''
    _BUF_SIZE = 1000
    
    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, options, params, lock):
        self.options = options
        self.params = params
        self.lock = lock
        self.reset_stats()
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def run(self, ibd_problems):
        '''Main call that outputs IBD segments of all hap pairs for all IbdProblem instances in the list
        ''ibd_problems''.
        
        Returns metadata: pair count and segment count for this object (since reset_stats() has
        been last called).'''
        self._reset_buffer()
        for ibd_problem in ibd_problems:
            self._pair_printout(ibd_problem, 'Started')
            # Do the core IBD calculation
            self.segments += im.ih.hap_segments(ibd_problem)
            self.buf += 1
            self.count_pairs += 1
            self._pair_printout(ibd_problem, 'Finished')
            if self.buf == PairProcessor._BUF_SIZE:
                # Buffer full, write results
                self._output_buffer()
                self._reset_buffer()
        if self.buf: self._output_buffer()  # Write remaining segment buffer, if exists
        return self.count_pairs, self.count_segments  # metadata

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def reset_stats(self):
        self.count_segments = 0
        self.count_pairs = 0

    def _reset_buffer(self):
        self.segments = im.segment.SegmentSet([])
        self.buf = 0

    def _output_buffer(self):
        '''Write results. Need a lock, otherwise output from the different processes
        is liable to getting all mixed up.'''
        self.lock.acquire()
        self.segments.save(sys.stdout)
        sys.stdout.flush()
        self.lock.release()
        self.count_segments += len(self.segments)

    def _pair_printout(self, ibd_problem, msg):
        if self.options.debug >= 2:    
            util.write_with_lock('Pair (%4d,%4d) (%4d,%4d): %s\n' % (ibd_problem.hap1[0], ibd_problem.hap1[1],
                                                                     ibd_problem.hap2[0], ibd_problem.hap2[1],
                                                                     msg), self.lock)

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s [flags] <phased-data-file> <kinship-file> <pair-list-file>\n\n' \
        'Locate IBD segments among a subset of sample in an NPZ phased data set.\n' \
        'Sample pairs are read from standard input. Segments are written to standard output.\n' \
        '\tphased-data-file - NPZ file containing phasing results\n' \
        '\tkinship-file - Sorted identity coefficient file\n' \
        '\tpair-list-file - Sorted identity coefficient file\n' \
        '\tout-file - File to output segments to\n' \
        '\n' \
        'Example:\n' \
        'phased-data-file = chr22/hutt.phased.npz\n' \
        'kinship-file = hutt.kinship\n' \
        'pair-list-file contains the lines\n' \
        '0 1\n' \
        '...\n' \
        '0 100\n' \
        '\n' \
        'Type ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-d', '--debug', type='int', dest='debug', default=0,
                      help='Debug Level (0=quiet; 1=summary; 2=full debug)')
    parser.add_option('-c', '--chunk-size', type='int', dest='chunk_size', default=10,
                      help='Number of pairs to pass to each process (batch processing within a process)')
    parser.add_option('-u', '--len-unit', type='str', dest='len_unit', default=im.PhaseParam().len_unit,
                      help='segment length to look for [cm|mbp]')
    parser.add_option('-l', '--min-len', type='str', dest='min_len', default=im.PhaseParam().min_len,
                      help='Minimum segment length to look for')
    parser.add_option('-p', '--processes', type='int', dest='num_processes', default=1,
                      help='Number of processes to spawn')
        
    options, args = parser.parse_args(sys.argv[1:])
    num_args = len(args)
    if num_args < 2 or num_args > 3:
        print usage
        sys.exit(1)
    phased_data_file, kinship_file = args[0:2]
    input_file = open(args[2], 'rb') if num_args >= 3 else sys.stdin
    
    try:
        start = time.time()
        
        # Create objects that are only needed once per all processes (problem is not passed in since it's
        # big and not serializable. Instead, we use it to create IbdProblems below).
        params = im.PhaseParam(kinship_file=kinship_file, debug=(options.debug >= 3), min_len=options.min_len)
        options.phased_data_file = phased_data_file
        problem = im.io.read_npz(options.phased_data_file)
        if options.debug >= 1:
            sys.stdout.write('Reading pairs from file %s, chunk_size %d\n' % (args[2], options.chunk_size))

        # Map phase: read sample pairs from file, create serializable IbdProblem objects for each,
        # and pass them to processes via a multiprocessing map
        manager = Manager()
        lock = manager.Lock()
        po = Pool(processes=options.num_processes)
        res = po.imap(process_pair, (([IbdProblem(problem, hap1, hap2, None, params) for hap1, hap2 in block],
                                      options, params, lock)
                                     for block in itemutil.grouper(read_pairs(input_file), options.chunk_size)))
        
        # Reduce phase
        total_pairs = sum(map(operator.itemgetter(0), res))
        total_segments = sum(map(operator.itemgetter(1), res))
        t = time.time() - start
        if options.debug >= 1:
            lock.acquire()
            sys.stdout.write('Total %d hap pairs, %d segments, %d processes\n' % (total_pairs, total_segments, options.num_processes))
            sys.stdout.write('Elapsed Time: %.3f sec (%.3f sec/pair)\n' % (t, t / total_pairs))
            lock.release()
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
