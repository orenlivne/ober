#!/usr/bin/env python
'''
============================================================
Functions related to the SNP LD graph, modules and frames. 

Created on November 21, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, util, networkx as nx, operator, itertools

#---------------------------------------------
# Methods
#---------------------------------------------
def load_ld_graph(file_name, threshold=1e-15):
    '''Load a chromosome's LD graph from an adjacency list text file.'''
    # print 'intersecting_genes', chrom, start, end
    g = nx.Graph()
    g.add_weighted_edges_from((i, j, w) for (i, j, w) in np.loadtxt(file_name,
                              dtype=[('i', 'S12'), ('j', 'S12'), ('r', 'f4')])
                              if w >= threshold)
    return g

def greedy_coloring(g, position):
    '''Color a graph for which the nodes correspond to a set of monotonically increasing positions on a line,
    and that approximately looks like a chain. position is a dictionary that maps nodes to positions.'''
    color = {}
    # Scan nodes by order of positions (left-to-right)
    #print '*' * 30
    for i in sorted(position, key=position.get):
        #print 'Scanning node', i, 'position', position[i]
        # Find an available color (does not appear in any of i's neighbors) of minimum value
        #print 'all neighbors', g.neighbors(i)
        #print color
        #print [color[j] for j in g.neighbors_iter(i) if color.has_key(j)]
        nbhr_colors = np.unique([color[j] for j in g.neighbors_iter(i) if color.has_key(j)])
        #print 'nbhr colors', [color[j] for j in g.neighbors_iter(i) if color.has_key(j)]
        #print 'Unique nbhr colors', nbhr_colors
        try:
            available_color = itertools.dropwhile(lambda (k, v): k == v, enumerate(nbhr_colors)).next()[1] - 1
            #print 'Found available color', available_color 
        except StopIteration:
            available_color = (nbhr_colors[len(nbhr_colors) - 1] + 1) if nbhr_colors.size else 0
            #print 'New color', available_color 
        color[i] = available_color
        #print 'assigning color', color[i], 'to', i 
        #print '*' * 30
    return color

def frames(chrom, snp_names, snp_dao, ld_dao, threshold=0.0):
    '''Split the list of SNPs 'snp_names' on chromosome 'chrom' into frames of independent SNPs. Note: if the original
    SNP set is sorted by base-pair locations, so will the output frames.'''
    if not isinstance(snp_names, np.ndarray):
        if not isinstance(snp_names, list): snp_names = list(snp_names)
        snp_names = np.array(snp_names)
    orig_indices = snp_names.argsort()
    ld_g = ld_dao.ld_graph(chrom, snps=snp_names, threshold=threshold)
    # SNPs that are independent of all other SNPs 
    independent = np.sort(orig_indices[np.searchsorted(snp_names[orig_indices], list(set(snp_names) - set(ld_g.nodes())))])
    return [np.sort(np.concatenate((independent, orig_indices[np.searchsorted(snp_names[orig_indices], list(frame))]))) 
            for frame in _compute_frames(ld_g, snp_dao)]
#     # Debugging: check that Illumina CytoSNP SNPs that are known to be in perfect LD are not in the same frame 
#     frames = [np.sort(np.concatenate((independent, orig_indices[np.searchsorted(snp_names[orig_indices], list(frame))]))) 
#             for frame in _compute_frames(ld_g, snp_dao)]
#     ind1 = np.where(snp_names == 'rs4819535')[0][0]
#     ind2 = np.where(snp_names == 'rs5748648')[0][0]
#     print [k for k, f in enumerate(frames) if ind1 in f]
#     print [k for k, f in enumerate(frames) if ind2 in f]
#     print map(len, frames)
#     return frames
        
def read_frames(in_file):
    '''Read frames from a text file.'''
    if isinstance(in_file, str):  # file_name is a file set prefix
        with open(in_file, 'rb') as f: return __read_frames(f)
    else: return __read_frames(in_file)  # f is a file handle
        
def write_frames(frames, out_file):
    '''Write frames to a text file.'''
    if isinstance(out_file, str):  # file_name is a file set prefix
        with open(out_file, 'wb') as f: __write_frames(frames, f) 
    else: __write_frames(frames, out_file)

####################################################################################
class Frames(object):
    def __init__(self, items):
        self._frames = util.mdict()
        for k, v in items: self._frames[k] = v

    #---------------------------------------------
    # Operators
    #---------------------------------------------
    def __repr__(self): return 'Frames[%s]' % (repr(self._frames),) 

    '''Equality of objects.'''
    def __eq__(self, other): return self.__key() == other.__key()
    def __ne__(self, other): return self.__key() != other.__key()

    def __key(self):
        '''Hash key.'''
        keys = sorted(self._frames.keys())
        return (keys, [[y.tolist() for y in self._frames[x]] for x in keys])

    '''Return the frame list of chromosome ''chrom''.'''
    def __getitem__(self, chrom): return self._frames[chrom]
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def iterkeys(self): return self._frames.iterkeys()
    def iteritems(self): return self._frames.iteritems()
    
    def to_array(self):
        keys = self._frames.keys()
        return np.array([keys, [self._frames[x] for x in keys]])

    @staticmethod
    def from_array(a): return Frames(zip(a[0], a[1]))

#---------------------------------------------
# Private Methods
#---------------------------------------------
'''Read frames from a text file (chromosome-to-frame-list dictionary)'''
__read_frames = lambda in_file: Frames((data[0], data[1:]) for data in  (np.fromstring(line, dtype=int, sep=' ') for line in in_file.readlines()))

def __write_frames(frames, out_file):
    '''Write frames to a text file.'''
    for chrom in frames.iterkeys():
        for f in frames[chrom]: out_file.write(repr(chrom) + ' ' + ' '.join('%d' % (x,) for x in f) + '\n')

def _compute_block_frames(g, snp_dao):
    '''Yield frame numbers (colors) within each SNP LD block using greedy coloring.'''
    # This is how to turn on database logging to debug query slowness:
#    logging.basicConfig()
#    logging.getLogger('sqlalchemy.engine').setLevel(logging.DEBUG)
#    logging.getLogger('sqlalchemy.pool.QueuePool').setLevel(logging.DEBUG)
    for block in nx.connected_components(g):
        snps = list(snp_dao.get_snps(block))
        position = dict((snp.name, snp.bp) for snp in snps)
#        print '#' * 50
#        print 'block', block 
#        print '#' * 50
#        print 'snps', snps 
#        print 'position', position 
        h = g.subgraph(block)
        c = greedy_coloring(h, position)
        frames = util.to_set_dict((v, k) for (k, v) in c.iteritems())
        yield frames

def _merge_block_frames(frames):
    '''Merge frames of different blocks to make them comparable in size.'''
    
    # Build a dictionary of (block #, frame # within block)- [uniquely identifying the block frame]
    # to-frame-size [SNPs]  
    block_frame_size = dict(sum(([((i, j), len(x)) for j, x in frame.iteritems()] for i, frame in enumerate(frames)), []))
    num_colors = max(len(frame.keys()) for frame in frames)

    # Loop over block frames in descending size order and append each to the currently smallest global
    # frame. This tends to equalize global frame sizes. 
    color_tot = np.zeros((num_colors,), dtype=np.uint)
    global_frame = [set() for _ in xrange(num_colors)]
    for frame, sz in sorted(block_frame_size.iteritems(), key=operator.itemgetter(1), reverse=True):
        color_of_min_size = np.argmin(color_tot)  # Could be accelerated by keeping color_tot sorted and updating it
        global_frame[color_of_min_size] |= frames[frame[0]][frame[1]]
        color_tot[color_of_min_size] += sz
    return global_frame

'''Yield frame numbers (colors) within each SNP LD block using greedy coloring. Frames are returned
in descending size order.'''
_compute_frames = lambda g, snp_dao: _merge_block_frames(list(_compute_block_frames(g, snp_dao)))
