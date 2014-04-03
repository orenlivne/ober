'''
============================================================
Naive node partitioning based on degrees (doesn't work very
well).

Created on March 26, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, itertools as it, networkx as nx, util

####################################################################################
class NaivePartitioner(object):
    '''Partitions a segment graph using a greedy approach (removes nodes with low degree from cliques;
    looks for extra edges from nearby segments first).'''    
    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, segments, bp_left, bp_right, min_degree, debug):
        self.segments = segments
        self.bp_left = bp_left
        self.bp_right = bp_right
        self.min_degree = min_degree
        self.debug = debug
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def partition(self, G):  # , initial_guess):
        return (c.nodes() for c in map(self._clean_clique, nx.connected_component_subgraphs(G)))
    
    def _clean_clique(self, c):
        '''Ideally, each IBD graph component is a clique (reflecting IBD transitivity).
        In practice, we have approximate transitivity, so attempt to complete the
        clique by looking for more IBD segment edges in the neighborhood of the SNP for
        nodes that are not connected to all others. If those nodes are still not connected
        to at least min_degree neighbors, remove them from the clique.
        If the clique is smaller than (min_degree+1), such nodes are always removed.'''
        max_degree = c.number_of_nodes() - 1
        nodes = c.nodes()
        orig_nodes = len(nodes) 
        unconnected = filter(lambda node: c.degree(node) < max_degree, nodes)
        if self.debug:
            self._writeln('\tComponent: nodes %d, edges %d, %s' % (c.number_of_nodes(), c.number_of_edges(), repr(util.occur_dict(c.degree().values()))))
        if unconnected:
            if self.debug:
                self._writeln('\t\tAll nodes %s' % (len(nodes),))
                self._writeln('\t\tUnconnected nodes %s' % (len(unconnected),))
            # Add extra edges
            extra_edges = [tuple(s.samples) for s in
                           it.chain.from_iterable(self.segments.find(self.bp_left, self.bp_right + 1, node)
                                                         for node in unconnected)
                           if s.samples[0] in nodes and s.samples[1] in nodes]
            # self._writeln('\t\textra_edges', extra_edges
            c.add_edges_from(extra_edges)
    #        if debug:
    #            self._writeln('\t\tnodes %d, edges %d, %s' % (c.number_of_nodes(), c.number_of_edges(), repr(util.occur_dict(c.degree().values())))
    
            # Delete still-connected nodes
            c_min_degree = min(max_degree, self.min_degree)
            unconnected = filter(lambda node: c.degree(node) < c_min_degree, nodes)
            if unconnected:
    #            if debug:
    #                self._writeln('\t\tRemoving nodes %s min_degree %d' % (repr(unconnected), c_min_degree)
                c.remove_nodes_from(unconnected)
            if self.debug:
                self._writeln('\t\tCleaned component: %d -> %d' % (orig_nodes, c.number_of_nodes()))
        return c
    
    def _writeln(self, s):    
        sys.stdout.write(s + '\n')
