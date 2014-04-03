#!/usr/bin/env python
'''
============================================================
Prune variants for LD. In each LD block, keep the variant
with the highest call rate. A sliding window implementation.

Created on January 21, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, linecache, networkx as nx, numpy as np, time
from optparse import OptionParser
from numpy.lib.function_base import corrcoef

def parse_command_line_args(argv):
    '''Parse and validate command-line arguments.'''
    prog = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <map-file> <genotype-file>\n\n' \
        'Prune variants for LD. In each LD block, keep the variant\n' \
        'with the highest call rate.\n' \
        'Input:\n' \
        'map-file: a file with a list of SNP names and call rates.\n' \
        'genotype-file: a file with genotypes in dosage format. One row per SNP. Matches map-file ordering. Missing genotypes should be encoded as -1, dosages as 0,1,2.\n\n' \
        'Output:\n' \
        'list of SNP IDs. Pruned with the specified r^2 threshold.\n' \
        'LD graph, saved in NPZ format to map-file.graph.npz.\n' \
        '\nType ''%s -h'' to display full help.' % (prog, prog)
    parser = OptionParser(usage=usage)
    parser.add_option('-s', '--separator',
                      default=' ',
                      help='Genotype-file delimiter [default: %default. Use the string ''\t'' for a tab]')
#    parser.add_option('-r', '--r2-threshold' , type='float', dest='r2_min', default=0.3,
#                      help='rsq threshold, continue calculating rsq until rsq is less than this value')
    parser.add_option('-p', '--pruning-threshold' , type='float', dest='r2_prune', default=0.99,
                      help='SNPs are partitioned into blocks such that rsq < this threshold for SNPs in different blocks.')
    parser.add_option('-v', '--debug'        , action='store_true'  , dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-f', '--print-frequency' , type='int', dest='frequency', default=100,
                      help='Print debugging information every this many SNPs')
    parser.add_option('-w', '--ld-window' , type='int', dest='ld_window', default=20,
                      help='At least these many neighbors of every will be considered in rsq calculations')
    parser.add_option('-k', '--ld-window_kb' , type='float', dest='ld_window_kb', default=500,
                      help='At least these many neighbors of every will be considered in rsq calculations')
    parser.add_option('-m', '--ld-window_shift' , type='int', dest='ld_window_shift', default=5,
                      help='Window shift step [#snps]')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 2:
        print usage
        sys.exit(1)
    # Support the tab delimiter - must convert it from the string '\t' to the character code \t
    # for the CSV reader to recognize it
    if options.separator == '\\t': options.separator = '\t'
    # Convert to kbp to bp
    options.ld_window_kb = int(options.ld_window_kb * 1000)
    return args, options

class GenotypeSet(object):
    '''Loads genotypes from a genotype file.'''

    def __init__(self, file_name, n):
        '''Initialize a genotype loader from the file \'file_name\'.'''
        self.file_name = file_name
#         self.n = n
#         
#        # Cache single-SNP moments
#         self.m1 = np.zeros((n,), dtype=float)
#         self.m2 = np.zeros((n,), dtype=float)
#         for i in xrange(n):
#             g = GenotypeSet.__parse_genotypes(linecache.getline(self.file_name, i + 1))
#             g = g[g >= 0] # Restrict to called genotypes 
#             self.m1[i] = sum(g)
#             self.m2[i] = sum(x * x for x in g)

    def __getitem__(self, i):
        '''Load genotype line i from the file f.'''
        return GenotypeSet.__parse_genotypes(linecache.getline(self.file_name, i + 1))

    @staticmethod  
    def __parse_genotypes(line):
        '''Parse a genotype line into a vector of genotype dosages.'''
        return np.array(map(float, line.split()))
    
def r2(x, y):
    '''Calculate the Pearson correlation coefficient between the vectors x and y using numpy.'''
    has_data = (x >= 0) & (y >= 0)
    return corrcoef(x[has_data], y[has_data])[0][1] ** 2

# def r2_mine(x, y): # Slower than numpy
#     '''Calculate the Pearson correlation coefficient between the vectors x and y.'''
#     has_data = (x >= 0) & (y >= 0)
#     x, y = x[has_data], y[has_data]
#     N = len(x)
#     mx, my = sum(x), sum(y)
#     sx, sy = sum(x * x), sum(y * y)
#     e = (N * sum(x * y) - mx * my)
#     return  (e * e) / ((N * sx - mx * mx) * (N * sy - my * my)) 

def r2_edgelist(file_name, n, options, snp, bp):    
    '''Yield the list of tuples (i,j,r^2(i,j)) in the file file_name, where i,j are line numbers.'''
    g = GenotypeSet(file_name, n)
    start_time = time.time()
    for i in xrange(0, n, options.ld_window_shift):
        Gi, j_right, bp_right = g[i], i + options.ld_window, bp[i] + options.ld_window_kb
        #if options.debug: print 'i', i, snp[i], bp[i], 'to', bp_right
        for j in xrange(i, n):
            rsq = r2(Gi, g[j])
            #if options.debug: print '\t', j, snp[j], rsq
            if j >= j_right or bp[j] >= bp_right: break  # If outside minimum window and small r^2, break
            yield i, j, rsq
        if options.debug:
            if i > 0 and i % options.frequency == 0:
                t = time.time() - start_time
                snps_per_sec = i/t
                sys.stderr.write('%d SNPs (%.2f%%) elapsed %s remaining %s (%.2e snp/min) %d neighbors)\n' % \
                                     (i, (100.*i)/n, 
                                      time.strftime('%H:%M:%S', time.gmtime(t)),
                                      time.strftime('%H:%M:%S', time.gmtime((n-i)/snps_per_sec)),
                                      60 * snps_per_sec, j - i))

def best_in_block(block, priority):
    '''Return the (a) best SNP (attains the highest priority) in a list block of SNPs.''' 
    return max((priority[snp], snp) for snp in block)[1]

'''Main program'''
if __name__ == '__main__':
    args, options = parse_command_line_args(sys.argv)
    map_file = args[0]
    #graph_file = map_file + '.graph.npz'
    
    # Map file format: snp bp call_rate is_func
    # Read SNP list into arrays; SNPs and priorities into a dictionary
    snp = [line.rstrip().split()[0] for line in open(map_file, 'rb')]
    bp = [int(line.rstrip().split()[1]) for line in open(map_file, 'rb')]
    # SNP's priority
    priority = dict((x[0], float(x[2])) for x in (line.rstrip().split() for line in open(map_file, 'rb')))
    # Is SNP functional (if so, the best in its block keeps a reference to its name so that we can find it later)
    is_func = dict((x[0], x[3] == '1') for x in (line.rstrip().split() for line in open(map_file, 'rb')))
    annotation = dict((x[0], x[4]) for x in (line.rstrip().split() for line in open(map_file, 'rb')))
    if options.debug: sys.stderr.write('Original #SNPs = %d\n' % (len(snp),))
    
    # Calculate r^2 and generate graph for the first time
    G = nx.Graph()
    G.add_nodes_from(snp)
    G.add_edges_from((snp[i], snp[j]) for i, j, rsq in r2_edgelist(args[1], len(snp), options, snp, bp) if rsq >= options.r2_prune)
    # Prune for LD, print best SNP in every LD block, save graph to file
    #if options.debug: sys.stderr.write('Saving LD graph to %s\n' % (graph_file,))

    c = nx.connected_components(G)
    if options.debug: sys.stderr.write('LD graph: %d nodes, %d edges, %d components\n' % (G.number_of_nodes(), G.number_of_edges(), len(c)))
    #np.savez(graph_file, G=np.array([G]))
    # For each block, print best snp in block (highest priority) and a comma-delimited list of functional
    # SNP in its block (format: name1:annotation1,name2:annotation2,...)
    for block in c:
        func_annotations=','.join(snp + ':' + annotation[snp] for snp in block if is_func[snp])
        print best_in_block(block, priority), func_annotations if func_annotations else '-'
