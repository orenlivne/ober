'''
============================================================
Print original graph weights and and affinities for edges
between sets of conflicting alleles (R1,R2)
identified at chr22:SNP#456 imputation.

Created on April 2, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, networkx as nx, itertools as it
from impute.dev.index_segments import conflict_cross_edges_init
# from scipy.sparse.csr import csr_matrix

G = conflict_cross_edges_init.G # singleton pattern
G = nx.connected_component_subgraphs(G)[0]
c, H, x, cliques = im.index.partition_amg.partition(G, full_output=True, theta=0.995)

cross_edges = lambda R1, R2: [(i, j, G[i][j]['weight'] if j in G[i] else 0.0, im.index.partition_amg.edge_affinity(x, [(G.nodes().index(i), G.nodes().index(j))])[0]) for (i, j) in it.product(R1, R2) if i < j]
print_cross_edges = lambda R1, R2: '\n'.join('(%4d,%d)  (%4d,%d)  %.3f  %.3f' % (a[0] + a[1] + (a[2], a[3])) for a in cross_edges(R1, R2))

R1 = [(633, 0), (567, 0), (807, 0), (257, 1), (199, 1), (571, 0), (807, 1), (406, 0), (1321, 1)]
R2 = [(789, 0)]

print 'R1 x R2'
print print_cross_edges(R1, R2)
print '-' * 50
print 'R1 x R1'
print print_cross_edges(R1, R1)
print '-' * 50
print 'R2 x R2'
print print_cross_edges(R2, R2)
