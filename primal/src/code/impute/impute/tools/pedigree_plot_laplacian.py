'''
============================================================
Home-grown Pedigree drawing tools. Graph Laplacian layout.
Not very successful so far.

Created on August 15, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, matplotlib.pyplot as plt, networkx as nx
import impute.tools.pedigree_tools as pt
from impute.data.Person import Person
from scipy.sparse.linalg.eigen.arpack.arpack import eigsh
from impute.tools.laplacian import laplacian

#---------------------------------------------
# Constants
#---------------------------------------------

#---------------------------------------------
# Methods
#---------------------------------------------
def plot(p):
    '''Draw pedigree.'''
    g = _marriage_graph(p.graph)
    pos = _layout_positions(p.graph, g)
    g = g.copy().to_undirected()
    
    # nodes
    node_size = 700
    alpha = 0.8
    
    data = [
            ([i for i in p.graph.nodes_iter() if p.person[i].node_type == Person.TYPE.NOT_GENOTYPED],
            'o', 'y'),
            ([i for i in p.graph.nodes_iter() if p.person[i].node_type == Person.TYPE.GENOTYPED],
            'o', 'r'),
            ([i for i in g if i < 0], '.', 'k')
            ]

    for (node_list, node_shape, node_color) in data:
        nx.draw_networkx_nodes(g, pos,
                               nodelist=node_list,
                               node_shape=node_shape,
                               node_color=node_color,
                               node_size=node_size, alpha=alpha)
    
    # edges
    nx.draw_networkx_edges(g, pos, width=1.0, alpha=0.5)
    # nx.draw_networkx_edges(G,pos,
    #                       edgelist=[(0,1),(1,2),(2,3),(3,0)],
    #                       width=8,alpha=0.5,edge_color='r')
    # nx.draw_networkx_edges(G,pos,
    #                       edgelist=[(4,5),(5,6),(6,7),(7,4)],
    #                       width=8,alpha=0.5,edge_color='m')
    
    # Some math labels
    labels = dict([(i, str(i)) for i in p.graph.nodes_iter()])
    nx.draw_networkx_labels(g, pos, labels, font_size=20)
    
    plt.axis('off')
    # plt.savefig('labels_and_colors.png') # save as png
    plt.show()  # display

#---------------------------------------------
# Private Methods
#---------------------------------------------

# Note: anything with a single underscore is a private method exposed to test classes. 

def _marriage_graph(g):
    '''Return an extended pedigree graph that includes marriage nodes and edges.
    Marriage node IDs are negative to avoid collision with real nodes.'''
    h = nx.DiGraph()
    h.add_nodes_from(g.nodes_iter())
    h.add_edges_from(g.edges_iter())
    marriage = -1  #  
    for (parents, children) in pt.families(g):
        h.add_edges_from((parent, marriage) for parent in parents)
        h.add_edges_from((marriage, child) for child in children)
        h.remove_edges_from((parent, child) for parent in parents for child in children)
        marriage -= 1 
    return h

def _layout_positions(g, g_extended=None):
    '''Return a dictionary that maps pedigree nodes to graph layout positions (x,y). Uses
    a layered spring model. Capable of handling both normal and extended pedigree graphs.'''
    g_extended = g_extended if g_extended else g
    # x-positions: spring model
    L = laplacian(g_extended, dtype=np.float)
    [_, v] = eigsh(L, 2, which='SM')
    x = v[:, 1]

    # Compute depths based on -original- graph
    y = pt.node_depth(g)
    for node in g_extended.nodes_iter():
        # Marriage node: 0.5 above child depth. Note that marriage nodes must have children
        if node < 0:
            y[node] = y[g_extended.successors_iter(node).next()] - 0.5
    
    ymax = max(d for d in y.itervalues())
    # Reverse the y-axis (generation 0 on top)
    return dict(zip(g_extended.nodes_iter(), zip(x, ((ymax - y[node]) for node in g_extended.nodes_iter()))))
