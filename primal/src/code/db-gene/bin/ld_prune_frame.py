#!/usr/bin/env python
'''
============================================================
Prune SNPs using LD R^2 data. Output the largest frame.

Created on August 8, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import itertools, sys, os, optparse, networkx as nx, numpy as np

#---------------------------------------------
# Private Methods
#---------------------------------------------
def to_set_dict(a):
    '''Convert a key-value-pair list a=[ (1, 'A'), (1, 'B'), (2, 'C') ] into a multi-valued dictionary,
    so that all values with the same key would be aggregated into a set:
    { 1: set(['A', 'B']), 2:(set['C'],) }. The dictionary is sorted by keys.'''
    return dict((key, set(v for (_, v) in pairs)) 
                for (key, pairs) in itertools.groupby(sorted(a), lambda pair: pair[0]))

def load_ld_graph(file_name, separator=' ', threshold=0, usecols=[0, 1, 2, 3, 4], skiprows=0):
    '''Create an LD graph from an adjacency list text file.'''
    # print 'intersecting_genes', chrom, start, end
    g = nx.Graph()
    edges = np.array([(i_bp, i, j_bp, j, w) for (i_bp, i, j_bp, j, w) in np.loadtxt(file_name, usecols=usecols, skiprows=skiprows, delimiter=separator,
                                                                    dtype=[('i_bp', 'i4'), ('i', 'S12'), ('j_bp', 'i4'), ('j', 'S12'), ('r', 'f4')])
                                                                    if w >= threshold],
                     dtype=[('i_bp', 'i4'), ('i', 'S12'), ('j_bp', 'i4'), ('j', 'S12'), ('r', 'f4')])
    position = dict(zip(edges['i'], edges['i_bp']))
    g.add_nodes_from([(snp, dict(position=bp)) for snp, bp in position.iteritems()])
    position = dict(zip(edges['j'], edges['j_bp']))
    g.add_nodes_from([(snp, dict(position=bp)) for snp, bp in position.iteritems()])
    g.add_weighted_edges_from(zip(edges['i'], edges['j'], edges['r']))
    return g

def largest_frame(ld_g):
    '''Split the list of SNPs 'snp_names' on chromosome 'chrom' into frames of independent SNPs. Note: if the original
    SNP set is sorted by base-pair locations, so will the output frames.'''
    snp_names = np.array(ld_g.nodes())
    orig_indices = snp_names.argsort()
    # SNPs that are independent of all other SNPs 
    independent = np.sort(orig_indices[np.searchsorted(snp_names[orig_indices], list(set(snp_names) - set(ld_g.nodes())))])
    frame = compute_frame(ld_g)
    return np.sort(np.concatenate((independent, orig_indices[np.searchsorted(snp_names[orig_indices], list(frame))]))) 

def greedy_coloring(g, position):
    '''Color a graph for which the nodes correspond to a set of monotonically increasing positions on a line,
    and that approximately looks like a chain. position is a dictionary that maps nodes to positions.'''
    color = {}
    # Scan nodes by order of positions (left-to-right)
    # print '*' * 30
    for i in sorted(position, key=position.get):
        # print 'Scanning node', i, 'position', position[i]
        # Find an available color (does not appear in any of i's neighbors) of minimum value
        # print 'all neighbors', g.neighbors(i)
        # print color
        # print [color[j] for j in g.neighbors_iter(i) if color.has_key(j)]
        nbhr_colors = np.unique([color[j] for j in g.neighbors_iter(i) if color.has_key(j)])
        # print 'nbhr colors', [color[j] for j in g.neighbors_iter(i) if color.has_key(j)]
        # print 'Unique nbhr colors', nbhr_colors
        try:
            available_color = itertools.dropwhile(lambda (k, v): k == v, enumerate(nbhr_colors)).next()[1] - 1
            # print 'Found available color', available_color 
        except StopIteration:
            available_color = (nbhr_colors[len(nbhr_colors) - 1] + 1) if nbhr_colors.size else 0
            # print 'New color', available_color 
        color[i] = available_color
        # print 'assigning color', color[i], 'to', i 
        # print '*' * 30
    # print 'Final colors', color 
    return color

def compute_block_frames(g):
    '''Yield frame numbers (colors) within each SNP LD block using greedy coloring.'''
    # This is how to turn on database logging to debug query slowness:
#    logging.basicConfig()
#    logging.getLogger('sqlalchemy.engine').setLevel(logging.DEBUG)
#    logging.getLogger('sqlalchemy.pool.QueuePool').setLevel(logging.DEBUG)
    for block in nx.connected_components(g):
        position = dict((snp, g.node[snp]['position']) for snp in block)
        # print '#' * 50
        # print 'block', l, block 
        # print '#' * 50
        # print 'position', [(x, position[x]) for x in sorted(position, key=position.get)]
        h = g.subgraph(block)
        c = greedy_coloring(h, position)
        frames = to_set_dict((v, k) for (k, v) in c.iteritems())
        yield frames

def merge_block_frames_to_largest(frames):
    '''Merge frames of different blocks. Return the largest frame.'''
    # Build a dictionary of (block #, frame # within block)- [uniquely identifying the block frame]
    # to-frame-size [SNPs]
    return set([x for f in frames for x in f[max([(len(v), k) for k, v in f.iteritems()])[1]]])

'''Yield frame numbers (colors) within each SNP LD block using greedy coloring. Frames are returned
in descending size order.'''
compute_frame = lambda g: merge_block_frames_to_largest(list(compute_block_frames(g)))

#---------------------------------------------
# Private Methods
#---------------------------------------------
def parse_command_line_args(argv):
    '''Parse and validate command-line arguments.'''
    PROGRAM = os.path.basename(argv[0])
    usage = 'Usage: %s [flags] \n\n' \
        'Read LD R^2 adjacency list from stdin, and write the largest frame of pruned SNPs\n' \
        'to write to stdout. Output the largest frame.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-s', '--separator', default=' ', dest='separator',
                      help='Delimiter [default: %default. Use the string ''\t'' for a tab]')
    parser.add_option('-n', '--skip-rows', type=int, default=0, dest='skiprows',
                      help='Number of header lines at the top of the file to skip')
    parser.add_option('-c', '--use-cols', default='0,1,2,3,4', dest='usecols',
                      help='Comma-separated list of three columns: bp_position_of_snp1, snp1, bp_position_of_snp2, snp2, R^2.')
    parser.add_option('-t', '--threshold', type='float', dest='threshold', default=0.0,
                      help='R^2 threshold above which (inclusive) SNPs are considered in LD.')
    options, args = parser.parse_args(argv[1:])
    if len(args) != 0:
        print usage
        sys.exit(1)
    # Support the tab delimiter - must convert it from the string '\t' to the character code \t
    # for the CSV reader to recognize it
    if options.separator == '\\t': options.separator = '\t'
    options.usecols = map(int, options.usecols.split(','))
    options.input = sys.stdin #if args[0] == 'stdin' else args[0] 
    options.output = sys.stdout #if args[1] == 'stdout' else args[1]
    return args, options

#---------------------------------------------
# Main Program
#---------------------------------------------           
if __name__ == '__main__':
    args, options = parse_command_line_args(sys.argv)
    # Load LD graph
    g = load_ld_graph(options.input, separator=options.separator, threshold=options.threshold,
                        usecols=options.usecols, skiprows=options.skiprows)
    # Calculate largest frame indices
    frame = largest_frame(g)
    # Translate back to snp names
    snp_names = np.array(g.nodes())
    snps = snp_names[frame]

    # Check that there are no edges within the frame
    for x in snps:
        for y in snps:
            if g.has_edge(x, y): raise ValueError('Frame did not pass quality control: edge exists between frame elements %s and %s' % (x, y))
    
    np.savetxt(options.output, snps, fmt='%s')
