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
import os, sys, impute as im, csv, optparse, traceback, util, networkx as nx, itertools as it, \
numpy as np, time, pickle, operator  # , cProfile, pstats
from impute.data.constants import ALLELES, MEGA_BASE_PAIR, PATERNAL, MATERNAL
from multiprocessing import Pool, Manager
from impute.ibd.index import partition_cast, partition_amg
from collections import OrderedDict

#---------------------------------------------
# Methods
#---------------------------------------------
def ibd_edge_weight(x, a, b, m0=2.0 * MEGA_BASE_PAIR, m_inf=0.4 * MEGA_BASE_PAIR, L_half_life=2.0 * MEGA_BASE_PAIR, mu=0.001):
    '''A rough estimate for the weight to be assigned to the graph edge between two IBD individuals at a SNP base-pair
    position x along a base-pair segment [a,b). M=margin half-life time [base-pairs].
    mu=confidence value at segment endpoints.
    The weight''s unit is ''effective mega-base-pairs'', and it is proportional to the segment length.'''
    # return mu + (1 - mu) * (b - a) * (1 - np.exp((a - x) / M)) * (1 - np.exp((x - b) / M)) / MEGA_BASE_PAIR
    M = m_inf + (m0 - m_inf) / (1 + (b - a) / L_half_life)
    return mu + (1 - mu) * (1 - np.exp((a - x) / M)) * (1 - np.exp((x - b) / M)) 

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
    if options.algorithm == 'amg': options.margin = 0.0
    if options.regions: options.regions = map(int, options.regions.split(','))
    options.save = (args[3] != '-')  # not options.snp_index is not None
    options.force_save_metadata = options.snp_index is not None or options.regions
    return args, options

def _save_index(region_info, result, save, out_dir):
    '''Save index to NPZ files, one per region. Regions are packed one after another in the result,
    so unpack them into separate blocks first.'''
    index = 0
    for r in region_info:
        num_snps = r['num_snps']
        block = result[index:index + num_snps]
        group_index = np.array(map(operator.itemgetter(0), block), np.int16)  # haplotype-to-IBD-group index
        groups = np.array([(np.array(group, dtype=np.int) if group is not None else None) for groups in map(operator.itemgetter(1), block) for group in groups])  # IBD-group-to-haplotype index
        if options.snp_index is not None: groups = groups[np.newaxis]
        # Only save if SNP index not specified
        if options.save:
            out_file = '%s/region-%d' % (out_dir, r['region'] * options.region_size)
            _writeln('Saving index to %s.npz' % (out_file,), None)
            np.savez(out_file, group_index=group_index, groups=groups)
        index += num_snps  # Advance pointer to start of next block in result list

'''Left and right margins around a BP position bp in a SNP BP array snp_np.'''
_bp_margins = lambda snp_bp, bp, margin: (max(snp_bp[0], bp - margin), min(snp_bp[-1] + 1, bp + margin))

def _bp_margin_index(snp_bp, start, value, direction):
    '''Look for the index enclosing a left/right margin value (direction=-1,1, repsectively)
    in a sorted SNP BP array.'''
    if direction < 0:
        while start >= 0 and value <= snp_bp[start]: start += direction
    else:
        while start < len(snp_bp) and value >= snp_bp[start]: start += direction
    return start - direction
    
def _estimated_call_rate(cliques, R, num_samples, num_haps):
    '''Estimate ideal genotype and allele call rates given the typed haplotype set R.'''
    imputable_haps = sum((c for c in cliques if set(c) & R), [])
    call_rate_hap = (100.*len(imputable_haps)) / num_haps
    call_rate_genotype = (100.*sum(1 for sample in xrange(num_samples) 
                                   if (sample, PATERNAL) in imputable_haps and
                                   (sample, MATERNAL) in imputable_haps)) / num_samples
    return call_rate_hap, call_rate_genotype

def _writeln(s, lock):
    '''Synchronize stdout line writing - main thread. Also works in single-processor mode.'''
    if lock: lock.acquire()
    sys.stdout.write(s + '\n')
    sys.stdout.flush()
    if lock: lock.release()

def _new_snp_processor(info, segments, options, lock):
    '''A factory method of SNP processors.'''
    processor = _SyncedRegionProcessor(info, segments, options, lock) if options.num_processes > 1 \
    else _SnpProcessor(info, segments, options)
    return processor

####################################################################################
class _SegmentCollection(object):
    '''A collection of chromosomal regions, each of which holds an IBD segment set.'''    
    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, info, segment_file, regions, options):
        '''Load segments from IBD segment list file.''' 
        # Do all I/O at once: load all segments into a list in memory
        if options.debug >= 1: self._writeln('Loading segments from file %s ...' % (segment_file,))
        segments = np.array(map(lambda line: map(int, line), csv.reader(open(segment_file, 'rb'), delimiter=' ', skipinitialspace=True)))
        # segments = np.loadtxt(open(segment_file, 'rb'), dtype=np.int) # slower

        # Create region metadata record array
        num_regions = len(regions)
        self.num_regions = num_regions
        r = np.zeros((num_regions,),
                     dtype=[('region', np.uint8),  # Region index
                            ('bp_start', np.uint),  # region start base-pair position
                            ('bp_stop', np.uint),  # region end base-pair position
                            ('bp_start_ext', np.uint),  # extended region (with margins) start base-pair position
                            ('bp_stop_ext', np.uint),  # extended region (with margins) stop base-pair position
                            ('snp_start', np.uint),  # region start base-pair SNP index
                            ('snp_stop', np.uint),  # region end base-pair SNP index
                            ('num_snps', np.uint)  # region end SNP index
                            ])
        snp_bp = info.snp['base_pair']
        for index, region in enumerate(regions):
            start_raw, stop_raw = region * options.region_size, min((region + 1) * options.region_size, info.num_snps)
            margin = int(options.margin * MEGA_BASE_PAIR)  # margin [bp]
            bp_left = _bp_margins(snp_bp, snp_bp[start_raw], margin)[0]
            bp_right = _bp_margins(snp_bp, snp_bp[stop_raw - 1] + 1, margin)[1]
            r[index] = (region, bp_left, bp_right,
                        _bp_margin_index(snp_bp, start_raw, bp_left, -1),
                        _bp_margin_index(snp_bp, stop_raw, bp_right, 1),
                        start_raw, stop_raw, stop_raw - start_raw)
        self.region_info = r
        
        # For each region, find intersecting segments and index them in a SmartSegmentSet object
        self._segments = OrderedDict() 
        min_len = options.min_len * MEGA_BASE_PAIR
        if options.debug >= 1: self._writeln('Total segments %d' % (len(segments),))
        for region_info in self.region_info:
            region, start, stop = region_info['region'], region_info['bp_start_ext'], region_info['bp_stop_ext']
            s = segments[(np.maximum(start, segments[:, 0]) < np.minimum(stop, segments[:, 1])) & 
                         (segments[:, 3] - segments[:, 2] >= min_len)]
            self._segments[region] = im.smart_segment_set.SmartSegmentSet.from_list(info.num_samples, s)
            if options.debug >= 1: self._writeln('Loaded region [%5d,%5d) %d' % (start, stop, self._segments[region].size))

    #---------------------------------------------
    # Operations
    #---------------------------------------------
    def __repr__(self):
        sizes = [r.size for r in self._segments.itervalues()]
        return 'SegmentCollection[regions=%d, total segments=%d, segments per region=%s]' % \
            (self.num_regions, sum(sizes), sizes)

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def _writeln(self, s):
        '''Synchronize stdout line writing.'''
        _writeln(s, None)
    
    def find(self, region, bp, sample):
        '''find all segments or overlapping the segment [start,stop] (units=base-pairs)
        between the haplotype ''source'' and haplotypes in ''target''. If target=None, searches in all
        haplotypes.'''
        return self._segments[region].find(bp, bp + 1, sample)

####################################################################################
class _SnpProcessor(object):
    '''Builds IBD clique index for SNPs.'''    
    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, info, segments, options):
        self.info = info
        self.snp_bp = info.snp['base_pair']
        self.segments = segments
        self.options = options
        # List of WGS haplotypes, for estimating call rates
        self.R = set(list(it.product(im.examples.wgs_sample_index(), ALLELES)))
        self.partitioner = self._new_partitioner(options)
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def process(self, (region, snp)):
        '''Build IBD index for SNPs snp from the segment collection segments.'''
        options = self.options
        start_time = time.time()
        if options.profile: self._writeln('\t[%-5d] start_time %s' % (snp, time.asctime()))
        info = self.info
        bp = self.snp_bp[snp]
        
        # Build IBD graph components from all segments intersecting the SNP's
        # base-pair position
        s_time = time.time()
        G = self._ibd_graph_from_segments(region, bp) 
        if options.profile: self._writeln('\t[%-5d] _ibd_graph_from_segments took %.2f s' % (snp, time.time() - s_time))
            
        if options.save and options.snp_index is not None:
            out_file = '%s/graph-%d.pickle' % (options.out_dir, options.snp_index)
            self._writeln('Saving graph to %s' % (out_file,))
            pickle.dump(G, open(out_file, 'wb'))
            # np.savez(out_file, G=np.array([G]))
    
        s_time = time.time()
        cliques, snp_group_index, snp_groups = self._index_graph(G) 
        if options.profile: self._writeln('\t[%-5d] index_graph took %.2f s' % (snp, time.time() - s_time))
        
        if options.debug >= 1:
            if options.profile: self._writeln('\t[%-5d] total %.2f s' % (snp, time.time() - start_time))
            # Estimate ideal genotype and allele call rates
            num_haps = 2 * info.num_samples
            num_clique_haps = sum(map(len, cliques))
            num_cliques = len(cliques)
            call_rate_hap, call_rate_genotype = _estimated_call_rate(cliques, self.R, info.num_samples, num_haps) 
            self._writeln('\tSNP %-5d bp %d: groups %d over %d haps; distinct haps %d call rate allele %.2f%% genotype %.2f%%' % \
                          (snp, bp, num_cliques, num_clique_haps, num_haps - num_clique_haps + num_cliques,
                           call_rate_hap, call_rate_genotype))
            if options.debug >= 2: self._writeln('')
            sys.stdout.flush()
    
        # Convert list-of-lists to an array. This assumes groups has at least one non-empty element.
        return snp_group_index, snp_groups

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def _new_partitioner(self, options):
        '''A functor for partitioning a graph.'''
        algorithm = options.algorithm
        if algorithm == 'cast': return partition_cast.partition
        elif algorithm == 'amg': return lambda G: partition_amg.partition(G, theta=options.threshold)
        else: raise ValueError('Unsupported clique partitioning algorithm ''%s''' % (algorithm,)) 

    def _index_graph(self, G):
        '''Convert an IBD graph into a union of cliques, encoded in the group_index and group output arrays.
        Clique (group) index is 1-based.''' 
        debug = self.options.debug
        components = nx.connected_component_subgraphs(G)
        if debug >= 2: print 'Original components', np.array([c.number_of_nodes() for c in components])
        # Partition each component into cliques by cleaning nodes that are insufficiently connected
        cliques = map(self.partitioner, components)
        if debug >= 2: print 'cliques', np.array(map(len, cliques))
        # Do not save single-node cliques - not useful for imputation.
        cliques = [y for x in cliques for y in x if len(y) >= 2]
    
        # Encode cliques in groups + group index
        group_index = np.zeros((self.info.num_samples, 2), dtype=np.int16)    
        # IBD groups are 1-based, add a dummy 0th-group
        groups = [None] * (len(cliques) + 1)
        for m, hap_list in enumerate(cliques, 1): 
            hap_array = np.array(hap_list)
            group_index[hap_array[:, 0], hap_array[:, 1]] = m
            groups[m] = hap_list
        return cliques, group_index, np.array(groups)
    
    def _ibd_graph_from_segments(self, region, bp):
        '''An internal call that forms a graph from IBD segments intersecting a base-pair position bp.'''
        segments, info = self.segments, self.info
        G = nx.Graph()
        G.add_weighted_edges_from((s.sample0, s.sample1, ibd_edge_weight(bp, s.bp_start, s.bp_stop))
                                   for s in it.chain.from_iterable(segments.find(region, bp, (sample, a)) 
                                                                   for sample in xrange(info.num_samples) for a in ALLELES))
        return G

    def _writeln(self, s):    
        sys.stdout.write(s + '\n')

####################################################################################
class _SyncedRegionProcessor(_SnpProcessor):
    '''Processes a SNP. For multi-processing.'''
    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, info, segments, options, lock):
        _SnpProcessor.__init__(self, info, segments, options)
        self.lock = lock

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def _writeln(self, s):
        '''Synchronize stdout line writing.'''
        _writeln(s, self.lock)     

class MyObject(object):
    def process(self, region, snp): return region, snp

####################################################################################
# Main program
####################################################################################
def __main(args, options):
    '''Main program - accepts an options struct.'''    
    # Parse and validate command-line arguments
    in_file, info_file, segment_file, out_dir = args
    options.out_dir = args[3]  # Useful shortcut

    try:
        # Initialize thread pool
        if options.num_processes > 1:
            manager = Manager()
            lock = manager.Lock()
        else: lock = None
        start = time.time()

        # Load SNP info
        info = im.io.read_info_npz(info_file)
        if options.debug >= 1:
            _writeln('haps %d, snps %d, region size %d snps, processes %d' % \
                      (2 * info.num_samples, info.num_snps, options.region_size, options.num_processes), lock)
        
        # Read list of regions to process from stdin/in_file. If empty, process all regions.
        # If a region index list is read, IT MUST BE CONTIGUOUS! (e.g. [3, 4, 5, 6])
        regions = map(int, ([options.snp_index / options.region_size] if options.snp_index is not None else
                            (options.regions if options.regions else
                            (sys.stdin if in_file == '-' else open(in_file, 'rb')).readlines())))
        num_regions = (info.num_snps + options.region_size - 1) / options.region_size
        if not regions: regions = range(num_regions)
        _writeln('regions ' + repr(regions) + ' num_regions ' + repr(num_regions) + 
                 ' segment threshold ' + repr(options.min_len) + ' Mbp algorithm ' + options.algorithm + ' margin ' + repr(options.margin), lock)
                
        # Save index metadata if first region is processed in this run
        if options.save: util.mkdir_if_not_exists(out_dir)
        if options.save and (options.force_save_metadata or 0 in regions):
            if options.debug >= 1: _writeln('Writing metadata to %s/metadata' % (out_dir,), lock)
            np.savez('%s/metadata' % (out_dir,), snp=info.snp, region_size=options.region_size)

#        segments = _SegmentCollection(info, segment_file, regions, options)
#        print segments

        # Map phase: process each SNP independently
        r = segments.region_info
#        snps = [(r['region'][0], options.snp_index - r['snp_start'][0])] if options.snp_index is not None else \
#        ((region, snp) for region, start_raw, stop_raw in zip(r['region'], r['snp_start'], r['snp_stop'])
#         for snp in xrange(start_raw, stop_raw))
        snp_processor = _new_snp_processor(info, segments, options, lock)

        if options.num_processes > 1:
            # Multi-process mode. SNPs are processed in parallel. 
            po = Pool(processes=options.num_processes)
            #a = MyObject()
            #result = po.map(process_snp, ((a, region, snp) for region, snp in snps))
            #result = po.map(process_snp, ((region, snp) for region, snp in snps))
            #result = po.map(snp_processor.process, ((region, snp) for region, snp in snps))
            result = po.map(snp_processor.process, ((region, snp) for region, snp in snps))
        else:
            # Single-process mode (sequential)
            #result = [snp_processor.process((region, snp)) for region, snp in snps] 
            for region in regions:
                result = [snp_processor.process(segment_file, region,((region, snp)) for region, snp in snps)] 
        print result
        
        # Reduce phase: organize results in array and save to npz files
        _save_index(segments.region_info, result, options.save, out_dir) 

        t = time.time() - start
        if options.debug >= 1: _writeln('Elapsed Time: %.3f sec (%.3f sec/region)' % (t, t / len(regions)), lock)
        if options.num_processes > 1: manager.shutdown()
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

#----------------- old code -----------------
# Load segments from segment file
#         name = 'index_segments'
#         cProfile.runctx('segments = _SegmentCollection(info, segment_file, regions, options); print segments',
#                         globals(), locals(), filename=name)
#         p = pstats.Stats(name).strip_dirs()
#         p.sort_stats('cumulative').print_stats(50)
