'''
============================================================
Print original graph weights and and affinities for edges
between sets of conflicting alleles (R1,R2)
identified at chr18:SNP#2609 imputation (found during
monogenics imputation at chr18:28659816).

Created on April 2, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, networkx as nx, itertools as it, pickle, os
from impute.data.constants import MEGA_BASE_PAIR
from impute.ibd.index.index_segments import ibd_edge_weight
from impute.ibd.index.index_segments import main
#import cProfile, pstats
# from impute.dev.index_segments import conflict_cross_edges_init
# from scipy.sparse.csr import csr_matrix

def print_some_weights():
    '''Debugging the clique [(46,1),(381,1),(655,1)].'''
    x = 28646066
    m0, m_inf, L_half_life = 2.0 * MEGA_BASE_PAIR, 0.4 * MEGA_BASE_PAIR, 2.0 * MEGA_BASE_PAIR
    
    a, b = 25705139, 78014581  # left (46, 1) (655, 1)
    M = m_inf + (m0 - m_inf) / (1 + (b - a) / L_half_life)
    print 'M = %f' % (M / MEGA_BASE_PAIR,)
    print ' x=%d [%d,%d] %.3f' % (x, a, b, ibd_edge_weight(x, a, b, m0=m0, m_inf=m_inf, L_half_life=L_half_life))
    a, b = 24393963, 28771842  # down
    M = m_inf + (m0 - m_inf) / (1 + (b - a) / L_half_life)
    print 'M = %f' % (M / MEGA_BASE_PAIR,)
    print ' x=%d [%d,%d] %.3f' % (x, a, b, ibd_edge_weight(x, a, b, m0=m0, m_inf=m_inf, L_half_life=L_half_life))
    a, b = 11385531, 78014581  # (46, 1) (381, 1) right
    M = m_inf + (m0 - m_inf) / (1 + (b - a) / L_half_life)
    print 'M = %f' % (M / MEGA_BASE_PAIR,)
    print ' x=%d [%d,%d] %.3f' % (x, a, b, ibd_edge_weight(x, a, b, m0=m0, m_inf=m_inf, L_half_life=L_half_life))
    
def load_graph(snp):
    print 'load_graph(%d)' % (snp,)
    return pickle.load(open(os.environ['OBER_OUT'] + '/requests/monogenic/ibd/graph-%d.pickle' % (snp,), 'rb'))

def index_segments():
    main('-',
         os.environ['OBER_OUT'] + '/phasing/chr18/hutt.phased.info.npz',
         os.environ['OBER_OUT'] + '/requests/monogenic/ibd/chr18/segments-2600-2700.out',
         '-',  # os.environ['OBER_OUT'] + '/requests/monogenic/ibd/chr18',
         debug=2, algorithm='amg', min_len=0.4, region_size=100, snp_index=2609)

#---------------------------------
# Main program
#---------------------------------
index_segments()

print_some_weights()
# G = conflict_cross_edges_init.G # singleton pattern
G = load_graph(2609)
G = nx.connected_component_subgraphs(G)[0]
c, H, x, cliques = im.index.partition_amg.partition(G, full_output=True, theta=0.995)

cross_edges = lambda R1, R2: [(i, j, G[i][j]['weight'] if j in G[i] else 0.0, im.index.partition_amg.edge_affinity(x, [(G.nodes().index(i), G.nodes().index(j))])[0]) for (i, j) in it.product(R1, R2) if i < j]
print_cross_edges = lambda R1, R2: '\n'.join('(%4d,%d)  (%4d,%d)  %.3f  %.3f' % (a[0] + a[1] + (a[2], a[3])) for a in cross_edges(R1, R2))

R1 = [(46, 1)]
R2 = [(972, 0), (1246, 1), (944, 0), (803, 1), (888, 0)]

print 'R1 x R2'
print print_cross_edges(R1, R2)
print '-' * 50
print 'R1 x R1'
print print_cross_edges(R1, R1)
print '-' * 50
print 'R2 x R2'
print print_cross_edges(R2, R2)
