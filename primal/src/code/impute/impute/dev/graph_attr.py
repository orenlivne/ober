#!/usr/bin/env python
'''
============================================================
A test of passing graph attributes by reference, not by
value. 

Created on March 13, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import networkx as nx

G = nx.DiGraph()
G.add_node(1, dict(a=[1, 2, 3]))
G.add_node(2, dict(a=None))
G.add_node(3, dict(a=[1, 2, 3]))

G.add_edges_from([(1, 2), (1, 3)])

G.node[2]['f'] = G.node[1]['a']
print G.node[1]
print G.node[2]

G.node[1]['a'][0] = 4
print G.node[1]
print G.node[2]

G.node[1]['a'] = [5, 2, 3]
print G.node[1]
print G.node[2]
