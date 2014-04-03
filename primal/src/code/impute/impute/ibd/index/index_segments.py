#!/usr/bin/env python
'''
============================================================
Convert an IBD segment file (produced by ibd_segments.py)
to a directory with NPZ files containing the segments and
IBD sets of each chromosomal region. Each region is indexed
by SNP (over the same SNP set used to generate the segments,
e.g., Affymetrix SNPs).

Created on March 15, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, impute as im, csv, optparse, traceback, util, networkx as nx, itertools as it, numpy as np, time, pickle, cProfile, pstats
from impute.data.constants import ALLELES, MEGA_BASE_PAIR, PATERNAL, MATERNAL
from multiprocessing import Pool, Manager
from impute.ibd.index import partition_cast, partition_amg
from impute.ibd.index.partition_naive import NaivePartitioner
from impute.data.problem import ProblemInfo

#---------------------------------------------
# Methods
#---------------------------------------------
def ibd_edge_weight(x, a, b, m0=2.0 * MEGA_BASE_PAIR, m_inf=0.4 * MEGA_BASE_PAIR, L_half_life=2.0 * MEGA_BASE_PAIR, mu=0.001):
    '''A rough estimate for the weight to be assigned to the graph edge between two IBD individuals at a SNP base-pair
    position x along a segment [a,b). M=margin half-life time [base-pairs]. mu=confidence value at segment endpoints.
    The weight''s unit is ''effective mega-base-pairs'', and it is proportional to the segment length.'''
    # return mu + (1 - mu) * (b - a) * (1 - np.exp((a - x) / M)) * (1 - np.exp((x - b) / M)) / MEGA_BASE_PAIR
    M = m_inf + (m0 - m_inf) / (1 + (b - a) / L_half_life)
    return mu + (1 - mu) * (1 - np.exp((a - x) / M)) * (1 - np.exp((x - b) / M)) 

def ibd_graph(info, segment_file, snp, min_len_mbp=0):
    '''Return the IBD segment graph at SNP index snp. Info is either a ProblemInfo object or file location.'''
    # Load segments intersecting the SNP's base-pair position
    if not isinstance(info, ProblemInfo):
        info = im.io.read_info_npz(info)  
    min_len = min_len_mbp * MEGA_BASE_PAIR
    segments = im.smart_segment_set.SmartSegmentSet.from_list(info.num_samples,
                                                    filter(lambda line: line[0] <= snp and snp < line[1] and line[3] - line[2] >= min_len,
                                                           map(lambda line: map(int, line), csv.reader(open(segment_file, 'rb'), delimiter=' ', skipinitialspace=True))))
    # Build IBD graph from segments 
    return _ibd_graph_from_segments(info, segments, info.snp['base_pair'][snp])
    
def _ibd_graph_from_segments(info, segments, bp):
    '''An internal call that forms a graph from IBD segments intersecting a base-pair position bp.'''
    G = nx.Graph()
#     for s in it.chain.from_iterable(segments.find(bp, bp + 1, (sample, a)) 
#                                     for sample in xrange(info.num_samples) for a in ALLELES):
#         x, a, b = bp, s.start, s.stop
#         print s.samples[0], s.samples[1],
#         print ' x=%d [%d,%d] %.3f' % (x, a, b, ibd_edge_weight(x, a, b))
    G.add_weighted_edges_from((s.sample0, s.sample1, ibd_edge_weight(bp, s.bp_start, s.bp_stop))
                              for s in it.chain.from_iterable(segments.find(bp, bp + 1, (sample, a)) 
                                                              for sample in xrange(info.num_samples) for a in ALLELES))
    return G

def index_ibd_graph(info, G, partitioner, debug=0):
    '''Convert an IBD graph into a union of cliques, encoded in the group_index and group output arrays.
    Clique (group) index is 1-based.''' 
    components = nx.connected_component_subgraphs(G)
    if debug >= 2:
        print 'Original components', np.array([c.number_of_nodes() for c in components])
    # Partition each component into cliques by cleaning nodes that are
    # insufficiently connected to other nodes.
    # TODO: use initial guess from previous SNP?
    # cliques = partitioner(G) #, initial_guess=cliques)
    cliques = map(partitioner, components)
    if debug >= 2:
        print 'cliques', np.array([np.array(map(len, c)) for c in cliques])
    # Do not save single-node cliques. They are not useful for imputation.
    cliques = [y for x in cliques for y in x if len(y) >= 2]

    # Encode cliques in groups + group index
    group_index = np.zeros((info.num_samples, 2), dtype=np.int16)    
    # IBD groups are 1-based, add a dummy 0th-group
    groups = [None] * (len(cliques) + 1)
    for m, hap_list in enumerate(cliques, 1): 
        hap_array = np.array(hap_list)
        group_index[hap_array[:, 0], hap_array[:, 1]] = m
        groups[m] = hap_list
    return cliques, group_index, groups
    
#---------------------------------------------
# Private Methods
#---------------------------------------------
def __parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s [flags] <input-file or -> <info-npz-file> <segment-file> <output-dir>\n\n' \
        'Convert an IBD segment file (produced by ibd_segments.py)\n' \
        'to a directory with NPZ files containing the segments and\n' \
        'IBD sets of each chromosomal region. Each region is indexed\n' \
        'by SNP (over the same SNP set used to generate the segments,\n' \
        'e.g., Affymetrix SNPs).\n' \
        '\n' \
        'Type ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-v', '--debug', type='int', dest='debug', default=0,
                      help='Print debugging messages (0: no messages 1: some messages 2: clique cleaning debugging 3: run only one test SNP')
    parser.add_option('-r', '--region-size', type='int', dest='region_size', default=200,
                      help='Number of SNPs in each region output file')
    parser.add_option('-w', '--regions', type='str', dest='regions', default=None,
                      help='Run these region numbers instead of reading from stdin or file (comma-separated list). If specified, metadata file is also written.')
    parser.add_option('-s', '--snp-index', type='int', dest='snp_index', default=None,
                      help='Run a single SNP index instead of reading from stdin or file. Overrides -w. If specified, metadata file is also written.')
    parser.add_option('-n', '--min-degree', type='int', dest='min_degree', default=4,
                      help='Minimum number of IBD-sharing samples that a sample must have to be in an IBD clique of size (min-degree+1) or larger')
    parser.add_option('-l', '--min-length', type='float', dest='min_len', default=2.0,
                      help='Minimum IBD segment length to consider [Mbp]')
    parser.add_option('-m', '--margin', type='float', dest='margin', default=0.8,
                      help='Base-pair margin around a SNP to look for neighboring IBD segments [Mbp]')
    parser.add_option('-p', '--processes', type='int', dest='num_processes', default=1,
                      help='Number of processes to spawn')
    parser.add_option('-a', '--algorithm', type='str', dest='algorithm', default='naive',
                      help='Quasi-clique partitioning algorithm: naive|amg.')
    parser.add_option('-t', '--threshold', type=float, dest='threshold', default=0.995,
                      help='Quasi-clique partitioning edge dropping affinity threshold.')
    parser.add_option('-f', '--profile', action='store_true', dest='profile', default=False,
                      help='Profile index building')

    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 4:
        sys.stdout.write(usage + '\n')
        sys.exit(1)
    return __validate_options(args, options)

def __validate_options(args, options):
    if options.algorithm == 'amg':
        options.margin = 0.0
    if options.regions:
        options.regions = map(int, options.regions.split(','))
    options.save = (args[3] != '-')  # not options.snp_index
    options.force_save_metadata = options.snp_index or options.regions
    return args, options

def _process_region((info, segment_file, out_dir, options, i, lock)):
    '''Process a region - this is what each processor runs.'''
    processor = SyncedRegionProcessor(lock) if options.num_processes > 1 else RegionProcessor()
    group_index, groups = processor.build_index(info, segment_file, options, i)
        
    # Only save if SNP index not specified
    if options.save:
        out_file = '%s/region-%d' % (out_dir, i * options.region_size)
        _writeln('Saving index to %s.npz' % (out_file,), lock)
        np.savez(out_file, group_index=group_index, groups=groups)
    return 0  # Dummy retval

def _process_region_profile((info, segment_file, out_dir, options, i, lock)):
    name = 'build_index_%d' % (options.region_size,)
    cProfile.runctx('_process_region((info, segment_file, out_dir, options, i, lock))',
                    globals(), locals(), filename=name)
    p = pstats.Stats(name).strip_dirs()
    p.sort_stats('cumulative').print_stats(50)

def _bp_margins(snp_bp, bp, margin):
    '''Return left and right margins around a BP position bp in a SNP BP array snp_np.'''
    return max(snp_bp[0], bp - margin), min(snp_bp[-1] + 1, bp + margin)

def _bp_margin_index(snp_bp, start, value, direction):
    '''Look for the index enclosing a left/right margin value (direction=-1,1, repsectively)
    in a sorted SNP BP array.'''
    if direction < 0:
        while start >= 0 and value <= snp_bp[start]:
            start += direction
    else:
        while start < len(snp_bp) and value >= snp_bp[start]:
            start += direction
    return start - direction

def _segment_attn(self, a, b, mu=0.9):
    '''Attenuate segment [a,b]''s size by a factor of mu if k is at the edge of the interval
    (no penalty in the middle), so that central segments are preferred.'''
    L = b - a 
    return L * (mu + (1 - mu) * np.abs(self.snp_bp - 0.5 * (a + b)) / (0.5 * L))

def _writeln(s, lock):
    '''Synchronize stdout line writing - main thread.'''
    if lock:     
        lock.acquire()
    sys.stdout.write(s + '\n')
    sys.stdout.flush()
    if lock:     
        lock.release()

####################################################################################
class RegionProcessor(object):
    '''Processes a chromosomal region. Not synchronized.'''    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def build_index(self, info, segment_file, options, i):
        '''Build IBD index for region #i of consecutive SNPs. Uses all segments that intersect
        the region + margin bp's on either side.'''
        # Calculate start and stop SNP indices of inner and augmented regions
        # (cf. inner-loop SNP margin)
        start_raw, stop_raw = i * options.region_size, min((i + 1) * options.region_size, info.num_snps)
        snp_bp = info.snp['base_pair']
        margin = int(options.margin * MEGA_BASE_PAIR)  # margin [bp]
        bp_left = _bp_margins(snp_bp, snp_bp[start_raw], margin)[0]
        bp_right = _bp_margins(snp_bp, snp_bp[stop_raw - 1] + 1, margin)[1]
        start = _bp_margin_index(snp_bp, start_raw, bp_left, -1)
        stop = _bp_margin_index(snp_bp, stop_raw, bp_right, 1)
        num_snps = stop_raw - start_raw 
                
        # Allocate haplotype-to-IBD-group index.
        # Assuming at most 2^16-1 distinct haplotypes exist in the data set. index=0 means not
        # in any IBD group. Note: IBD group indices are 1-based.
        group_index = np.zeros((num_snps, info.num_samples, 2), dtype=np.int16)    
        # Allocate IBD-group-to-haplotype index (list of lists for each SNP).
        groups = [list() for _ in xrange(num_snps)]

        # Load all segments intersecting the augmented region [start,stop)
        min_len = options.min_len * MEGA_BASE_PAIR
        if options.debug >= 1:
            self._writeln('Loading region [%d,%d) augmented [%d,%d) ...' % (start_raw, stop_raw, start, stop))
        segments = im.smart_segment_set.SmartSegmentSet.from_list(info.num_samples,
                                                        filter(lambda line: max(start, line[0]) < min(stop, line[1]) and line[3] - line[2] >= min_len,
                                                               map(lambda line: map(int, line), csv.reader(open(segment_file, 'rb'), delimiter=' ', skipinitialspace=True))))
        if options.debug >= 1:
            self._writeln('Loaded region [%d,%d) augmented [%d,%d) segments %d' % (start_raw, stop_raw, start, stop, segments.size))
            num_haps = 2 * info.num_samples
            R = set(list(it.product(im.examples.wgs_sample_index(), ALLELES)))
            
        partitioner = self._new_partitioner(options.algorithm, segments, bp_left, bp_right, options)
        # cliques = None
        
        # Loop over SNPs in region and append IBD groups to index
        for k, snp in enumerate([options.snp_index], options.snp_index - start_raw) \
        if options.snp_index else enumerate(xrange(start_raw, stop_raw)):
            start_time = time.time()
            if options.profile: self._writeln('\t[%-5d] start_time %s' % (snp, time.asctime()))

            bp = snp_bp[snp]  # SNP's base-pair position
            # Left and right margins to look for extra IBD segments around the SNP
            bp_left, bp_right = _bp_margins(snp_bp, bp, margin)

            # Build IBD graph components from all segments intersecting the SNP's
            # base-pair position
            s_time = time.time()
            G = _ibd_graph_from_segments(info, segments, bp) 
            # G = nx.from_edgelist((tuple(s.samples) for s in it.chain.from_iterable(segments.find(bp, bp + 1, (sample, a)) 
            #                                                                              for sample in xrange(info.num_samples) for a in ALLELES)))            
            if options.profile: self._writeln('\t[%-5d] _ibd_graph_from_segments took %.2f s' % (snp, time.time() - s_time))
                
            if options.save and options.snp_index is not None:
                out_file = '%s/graph-%d.pickle' % (options.out_dir, options.snp_index)
                self._writeln('Saving graph to %s' % (out_file,))
                pickle.dump(G, open(out_file, 'wb'))
                # np.savez(out_file, G=np.array([G]))
                
            s_time = time.time()
            cliques, snp_group_index, snp_groups = index_ibd_graph(info, G, partitioner, debug=options.debug) 
            if options.profile: self._writeln('\t[%-5d] index_ibd_graph took %.2f s' % (snp, time.time() - s_time))
            group_index[k] = snp_group_index
            groups[k] = snp_groups
            
            if options.debug >= 1:
                num_clique_haps = sum(map(len, cliques))
                num_cliques = len(cliques)
                # Estimate ideal genotype and allele call rates
                imputable_haps = sum((c for c in cliques if set(c) & R), [])
                call_rate_hap = (100.*len(imputable_haps)) / num_haps
                call_rate_genotype = (100.*sum(1 for sample in xrange(info.num_samples) 
                                               if (sample, PATERNAL) in imputable_haps and
                                               (sample, MATERNAL) in imputable_haps)) / info.num_samples
                if options.profile: self._writeln('\t[%-5d] total %.2f s' % (snp, time.time() - start_time))
                self._writeln('\tSNP %-5d bp %d (%d--%d): groups %d over %d haps; distinct haps %d call rate allele %.2f%% genotype %.2f%%' % \
                              (snp, bp, bp_left, bp_right, num_cliques, num_clique_haps, num_haps - num_clique_haps + num_cliques,
                               call_rate_hap, call_rate_genotype))
                if options.debug >= 2:
                    self._writeln('')
                sys.stdout.flush()


        # Convert list-of-lists to an array. This assumes groups has at least one non-empty element.
        return group_index, np.array([[np.array(x) for x in group] for group in groups])
                
    def _new_partitioner(self, algorithm, segments, bp_left, bp_right, options):
        '''A functor for partitioning a graph.'''
        if algorithm == 'naive': 
            p = NaivePartitioner(segments, bp_left, bp_right, options.min_degree, (options.debug >= 2))
            return p.partition
        elif algorithm == 'cast':
            return partition_cast.partition
        elif algorithm == 'amg':
            return lambda G: partition_amg.partition(G, theta=options.threshold)
        else:
            raise ValueError('Unsupported clique partitioning algorithm ''%s''' % (algorithm,)) 

    def _writeln(self, s):    
        sys.stdout.write(s + '\n')
                    
####################################################################################
class SyncedRegionProcessor(RegionProcessor):
    '''Processes a chromosomal region. For multi-processing.'''
    
    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, lock):
        RegionProcessor.__init__(self)
        self.lock = lock

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def _writeln(self, s):
        '''Synchronize stdout line writing.'''
        _writeln(s, self.lock)     

####################################################################################
# Main program
####################################################################################
def __main(args, options):
    '''Main program - accepts an options struct.'''    
    # Parse and validate command-line arguments
    in_file, info_file, segment_file, out_dir = args
    options.out_dir = args[3]  # Useful shortcut

    try:
        if options.num_processes > 1:
            manager = Manager()
            lock = manager.Lock()
        else:
            lock = None
        start = time.time()

        # Load SNP info
        info = im.io.read_info_npz(info_file)
        if options.debug >= 1:
            _writeln('haps %d, snps %d, region size %d snps, processes %d' % \
                      (2 * info.num_samples, info.num_snps, options.region_size, options.num_processes), lock)
        
        # Read list of regions to process from stdin/in_file. If empty, process all regions
        regions = map(int, ([options.snp_index / options.region_size] if options.snp_index else
                            (options.regions if options.regions else
                            (sys.stdin if in_file == '-' else open(in_file, 'rb')).readlines())))
        num_regions = (info.num_snps + options.region_size - 1) / options.region_size
        if not regions:
            regions = range(num_regions)
        _writeln('regions ' + repr(regions) + ' num_regions ' + repr(num_regions) + 
                  ' segment threshold ' + repr(options.min_len) + ' Mbp algorithm ' + options.algorithm + ' margin ' + repr(options.margin), lock)
        
        # Process each SNP region [start,stop) independently
        if options.save:
            util.mkdir_if_not_exists(out_dir)
        
        # Save index metadata, if processing the first region
        if options.save and (options.force_save_metadata or 0 in regions):
            if options.debug >= 1:
                _writeln('Writing metadata to %s/metadata' % (out_dir,), lock)
            np.savez('%s/metadata' % (out_dir,), snp=info.snp, region_size=options.region_size)

        process = _process_region_profile if options.debug >= 2 else _process_region
        if options.num_processes > 1:
            # Multi-process mode. Map phase:build and save regional index files 
            po = Pool(processes=options.num_processes)
            po.map(process, ((info, segment_file, out_dir, options, i, lock) for i in (i for i in regions if i >= 0 and i < num_regions)))
        else:
            # Single-process mode.
            for i in (i for i in regions if i >= 0 and i < num_regions):
                process((info, segment_file, out_dir, options, i, None))
            
        # Reduce phase - nothing to do here
        t = time.time() - start
        if options.debug >= 1:
            _writeln('Elapsed Time: %.3f sec (%.3f sec/region)' % (t, t / len(regions)), lock)
        if options.num_processes > 1:
            manager.shutdown()
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)

def main(in_file, info, segment_file, out_dir, **kwargs):
    '''Main program - reads options from a method parameter dictionary.'''
    # Default options
    options = util.Struct(num_processes=1, region_size=10, regions=None, snp_index=None, profile=False,
                          min_degree=4, min_len=2.0, margin=0.8, debug=0, algorithm='naive', threshold=0.995)
    # Override with passed arguments
    options.update(**kwargs)
    
    args = in_file, info, segment_file, out_dir
    args, options = __validate_options(args, options)
    return __main(args, options)
    
if __name__ == '__main__':
    '''Main program - accepts CLI arguments.'''
    args, options = __parse_command_line_args()
    __main(args, options)
