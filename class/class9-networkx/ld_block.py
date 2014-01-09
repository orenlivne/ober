#!/usr/bin/env python
'''
============================================================
Group SNPs into LD blocks. 
Convert an LD edge list to LD blocks. stdin contains a list 
of SNP pairs assumed to be in LD. The program lists the SNPs
and their 1-based LD block number (=connected component
number in the LD graph) into stdout. Block ordering is 
arbitrary.

Created on July 9, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, networkx as nx, sys

'''Load a chromosome's LD graph from an adjacency list text file. Pairs are assumed to be "in LD".
Note: in general, a weight (r^2) could and should be assigned to an edge.'''
read_graph = lambda file_name: nx.from_edgelist(((i, j) for (i, j) in np.loadtxt(file_name, dtype=[('i', 'S12'), ('j', 'S12')])), nx.Graph())

'''Write a list of SNP LD blocks to file (one line per SNP: SNP name, 1-based block number).'''
write_blocks = lambda blocks, out: np.savetxt(out, [(i, k) for k, block in enumerate(blocks, 1) for i in block], fmt='%s')

'''Main program'''
if __name__ == '__main__':
    write_blocks(nx.connected_components(read_graph(sys.stdin)), sys.stdout)
