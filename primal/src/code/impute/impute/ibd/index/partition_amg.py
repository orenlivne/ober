'''
============================================================
Partition a weighted undirected graph into quasi-cliques
using a two-level Algebraic Multigrid (AMG) coarsening
(based on algebraic distances). 

@see Clustering Gene Expression Patterns,
Amir Ben-Dor, Ron Shamir, Zohar Yakhiniz, 1999

Created on March 26, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, itertools, operator, networkx as nx
from impute.tools.laplacian import adjacency, diag

#---------------------------------------------
# Methods
#---------------------------------------------
def partition(G, K=5, nu=5, omega=0.7, theta=0.9, theta_weight=0.85, full_output=False):
    '''Partition the weighted undirected graph G into quasi-cliques using AMG coarsening.
    
    theta=affinity threshold. theta_weight=weight threshold. An edge must have
    affinity > theta and weight > theta_weight to stay in the graph.'''
    # Label nodes as a sequence integers. Otherwise A's rows ordering seems to be inconsistent
    # with G's node ordering
    original_nodes = G.nodes()
    G = nx.convert_node_labels_to_integers(G)
    
    # Original edge weights
    A = adjacency(G)
    ij = edge_list(A)
    w = A.data
    V = G.nodes()

    if theta > 0:
        # Calculate affinities
        x = test_vectors(A, K, nu, omega=omega)
        c = edge_affinity(x, ij)
        
        # New: Remove edges with small affinities from graph, unless they have a large edge weight w and
        # not a significantly smaller affinity. Fixes some imputation discordances in monogenic variant validation  
        # k = np.where((c < theta) & ((w < theta) | ((w >= theta) & (c < theta2))))[0]
        k = np.where((c < theta) | (w < theta_weight))[0]
    else:
        k = np.where(w < theta_weight)[0]
    
    # Old: Remove edges with small affinities from graph
    # k = np.where(c < theta)[0]
    # print len(ij), len(k), len(np.where(c < theta)[0])
    
    G2 = G.copy()
    e = [(V[ij[t][0]], V[ij[t][1]]) for t in k]
    G2.remove_edges_from(e)
#    print x
#    print K, nu, omega, theta
#     print 'V', V
#     print 'e', e
#     print len(e)
#     print 'c[0,39]', edge_affinity(x, [(0, 39)])
#     print 'c[e]', edge_affinity(x, e)
#     print G.number_of_nodes(), G.number_of_edges()
#     print G2.number_of_nodes(), G2.number_of_edges()
    
    # Convert back to original coordinates
    cliques = [map(original_nodes.__getitem__, clique) for clique in nx.connected_components(G2)]
#    print cliques
    if full_output:
        return c, G2, x, cliques 
    else:
        return cliques

def edge_list(A, rows=None):
    '''i-j edge list of selected rows or all rows (if None) of a CSR sparse matrix A.''' 
    return [(i, j) for i in (rows if rows is not None else xrange(A.shape[0])) for j in A.indices[A.indptr[i]:A.indptr[i + 1]]]

def edge_affinity(x, ij_list):
    return np.array([affinity_similarity(x[i], x[j]) for i, j in ij_list])

def affinity_distance(x, y):
    '''Return the affinity distance between test vector nodal value arrays x and y.'''
    return 1 - sum(x * y) ** 2 / (sum(x * x) * sum(y * y))

def affinity_similarity(x, y):
    '''Return the affinity similarity between test vector nodal value arrays x and y.'''
    return sum(x * y) ** 2 / (sum(x * x) * sum(y * y))

def to_logscale(d):
    '''Scale affinity similarity to a log scale similarity.''' 
    x = -np.log(1 - d)
    return x / (1 + x)

def test_vectors(A, K, nu, omega=0.7):
    '''Return K test vectors for a graph adjacency matrix A, obtained by nu omega-Jacobi
    relaxation sweeps, starting from random[-1,1].'''
    n = A.shape[0]
    # Initial test vectors
    np.random.seed(0)  # For reproducibility
    x = 2 * np.random.random((n, K)) - 1
    # Perform nu omega-Jacobi relaxations
    return relax(A, x, nu, omega=omega)

def relax(A, x, nu, omega=0.7):
    '''Perform nu omega-Jacobi relaxations on Ax=0, starting from the x array passed in,
    and returning the results.'''
    omega_d_inv = omega * np.tile(1. / diag(A), (x.shape[1],))
    for _ in xrange(nu):
        x = (1 - omega) * x + omega_d_inv * (A * x)
    return x

def cast_clustering(A, t=0.5, max_iter=1000):
    '''Implementation of the CAST quasi-clique partitioning algorithm.'''
    n = A.shape[0]
    C = []  # The collection of closed clusters
    U = set(range(n))  # Elements not yet assigned to any cluster
    count = 0
    while U and count < max_iter:
        count += 1
        Copen = []  # Start a new cluster
        a = dict((x, 0) for x in U)  # Reset affinity

        while True:  # While U contains a high affinity element, execute ADD step
            if not U:
                break
            u, au = max(((k, v) for k, v in a.iteritems() if k in U), key=operator.itemgetter(1))
            print 'ADD', u, au, t, len(Copen), '|U|', len(U)
            if au < t * len(Copen):
                break
            Copen.append(u)
            U.remove(u)
            for x in itertools.chain(U, Copen):
                a[x] = a[x] + A[x, u]
        
        while True:  # While Copen contains a low affinity element, execute REMOVE step
            if not Copen:
                break
            u, au = max(((k, v) for k, v in a.iteritems() if k in Copen), key=operator.itemgetter(1))
            print 'REMOVE', u, au, t, len(Copen)
            if au >= t * len(Copen):
                break
            Copen.remove(u)
            U.add(u)
            for x in itertools.chain(U, Copen):
                a[x] = a[x] - A[x, u]
        
        # Close the cluster
        C.append(Copen) 
    return C
