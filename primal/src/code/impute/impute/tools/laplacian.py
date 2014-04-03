'''
============================================================
Graph Laplacian functions.

Created on August 15, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, networkx as nx
from scipy.sparse.construct import spdiags

#---------------------------------------------
# Constants
#---------------------------------------------

#---------------------------------------------
# Methods
#---------------------------------------------
def adjacency(g, nodelist=None, weight='weight', dtype=None, form='csr'):
    '''Return the sparse adjacency matrix corresponding to a graph g. If g is undirected, returns
    a symmetric matrix.'''
    return nx.to_scipy_sparse_matrix(g, nodelist=nodelist, weight=weight, dtype=dtype, format=form)

def sym_adjacency(g, nodelist=None, weight='weight', dtype=None, form='csr'):
    '''Return the symmetrized sparse adjacency matrix corresponding to a graph g.'''
    A = nx.to_scipy_sparse_matrix(g, nodelist=nodelist, weight=weight, dtype=dtype, format=form)
    A = A + A.T
    return A

def diag(A, form='csr'):
    '''Return the diagonal of the Laplacian matrix corresponding to the adjacency matrix A.'''
    return np.asarray(A.sum(axis=1))

def laplacian(g, nodelist=None, weight='weight', dtype=None, form='csr'):
    '''Return the sparse graph Laplacian of a graph g.'''
    A = sym_adjacency(g, nodelist=nodelist, weight=weight, dtype=dtype, form=form)
    n = g.number_of_nodes()
    D = spdiags(np.asarray(A.sum(axis=1).T), 0, n, n, format=form)
    return D - A
