#!/usr/bin/env python
'''
============================================================
Draw a pedigree. This is a python translation of the
pedfiddler v0.6 algorithm by Loredo-Osti and Morgan
(http://www.stat.washington.edu/thompson/Genepi/Pedfiddler.shtml)
See license therein.

This code contains reproduction (cut-and-paste) of relevant
code from the NetworkX project (http://networkx.lanl.gov/index.html)
See license therein.

Input: a text file whose row format is
<id> <father_id> <mother_id> <sex> <affection_status> <label>

ID coding: missing=0, regular ID=any positive number
Sex coding: unknown=0, male=1, female=2
Affection status coding: unknown=0, normal=1, affected=2

Output: encapsulated postscript (EPS) file

NOTE: At this time, we do not support polygamous pedigrees.

Created on August 20, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
__author__ = 'Oren Livne <livne@uchicago.edu>'
__version__ = '1.0.0'

import sys, os, csv, time
from optparse import OptionParser, Option, OptionValueError
from fractions import gcd
from copy import deepcopy
from collections import defaultdict

########################################################################
# Constants
########################################################################
'''Index of element not found in a list.''' 
NOT_FOUND = -1

'''Indicates missing data in adjacency lists.'''
MISSING = 0

'''Enumerated types representing our coding conventions'''
class NODE_TYPE: DUMMY, PERSON, MARRIAGE = range(3)
class SEX: UNKNOWN, MALE, FEMALE = range(3)
class STATUS: UNKNOWN, NORMAL, AFFECTED = range(3)

'''
============================================================
Graph, from networkx
============================================================
'''
"""Base class for undirected graphs.

The Graph class allows any hashable object as a node
and can associate key/value attribute pairs with each undirected edge.

Self-loops are allowed but multiple edges are not (see MultiGraph).

For directed graphs see DiGraph and MultiDiGraph.
"""
#    Copyright (C) 2004-2011 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

__author__ = """\n""".join(['Aric Hagberg (hagberg@lanl.gov)',
                            'Pieter Swart (swart@lanl.gov)',
                            'Dan Schult(dschult@colgate.edu)'])

class Graph(object):
    """
    Base class for undirected graphs.

    A Graph stores nodes and edges with optional data, or attributes.

    Graphs hold undirected edges.  Self loops are allowed but multiple
    (parallel) edges are not.

    Nodes can be arbitrary (hashable) Python objects with optional
    key/value attributes.

    Edges are represented as links between nodes with optional
    key/value attributes.

    Parameters
    ----------
    data : input graph
        Data to initialize graph.  If data=None (default) an empty
        graph is created.  The data can be an edge list, or any
        NetworkX graph object.  If the corresponding optional Python
        packages are installed the data can also be a NumPy matrix
        or 2d ndarray, a SciPy sparse matrix, or a PyGraphviz graph.
    attr : keyword arguments, optional (default= no attributes)
        Attributes to add to graph as key=value pairs.

    See Also
    --------
    DiGraph
    MultiGraph
    MultiDiGraph

    Examples
    --------
    Create an empty graph structure (a "null graph") with no nodes and
    no edges.

    >>> G = Graph()

    G can be grown in several ways.

    **Nodes:**

    Add one node at a time:

    >>> G.add_node(1)

    Add the nodes from any container (a list, dict, set or
    even the lines from a file or the nodes from another graph).

    >>> G.add_nodes_from([2,3])
    >>> G.add_nodes_from(range(100,110))
    >>> H=Graph()
    >>> H.add_path([0,1,2,3,4,5,6,7,8,9])
    >>> G.add_nodes_from(H)

    In addition to strings and integers any hashable Python object
    (except None) can represent a node, e.g. a customized node object,
    or even another Graph.

    >>> G.add_node(H)

    **Edges:**

    G can also be grown by adding edges.

    Add one edge,

    >>> G.add_edge(1, 2)

    a list of edges,

    >>> G.add_edges_from([(1,2),(1,3)])

    or a collection of edges,

    >>> G.add_edges_from(H.edges())

    If some edges connect nodes not yet in the graph, the nodes
    are added automatically.  There are no errors when adding
    nodes or edges that already exist.

    **Attributes:**

    Each graph, node, and edge can hold key/value attribute pairs
    in an associated attribute dictionary (the keys must be hashable).
    By default these are empty, but can be added or changed using
    add_edge, add_node or direct manipulation of the attribute
    dictionaries named graph, node and edge respectively.

    >>> G = Graph(day="Friday")
    >>> G.graph
    {'day': 'Friday'}

    Add node attributes using add_node(), add_nodes_from() or G.node

    >>> G.add_node(1, time='5pm')
    >>> G.add_nodes_from([3], time='2pm')
    >>> G.node[1]
    {'time': '5pm'}
    >>> G.node[1]['room'] = 714
    >>> G.nodes(data=True)
    [(1, {'room': 714, 'time': '5pm'}), (3, {'time': '2pm'})]

    Warning: adding a node to G.node does not add it to the graph.

    Add edge attributes using add_edge(), add_edges_from(), subscript
    notation, or G.edge.

    >>> G.add_edge(1, 2, weight=4.7 )
    >>> G.add_edges_from([(3,4),(4,5)], color='red')
    >>> G.add_edges_from([(1,2,{'color':'blue'}), (2,3,{'weight':8})])
    >>> G[1][2]['weight'] = 4.7
    >>> G.edge[1][2]['weight'] = 4

    **Shortcuts:**

    Many common graph features allow python syntax to speed reporting.

    >>> 1 in G     # check if node in graph
    True
    >>> [n for n in G if n<3]   # iterate through nodes
    [1, 2]
    >>> len(G)  # number of nodes in graph
    5
    >>> G[1] # adjacency dict keyed by neighbor to edge attributes
    ...            # Note: you should not change this dict manually!
    {2: {'color': 'blue', 'weight': 4}}

    The fastest way to traverse all edges of a graph is via
    adjacency_iter(), but the edges() method is often more convenient.

    >>> for n,nbrsdict in G.adjacency_iter():
    ...     for nbr,eattr in nbrsdict.items():
    ...        if 'weight' in eattr:
    ...            (n,nbr,eattr['weight'])
    (1, 2, 4)
    (2, 1, 4)
    (2, 3, 8)
    (3, 2, 8)
    >>> [ (u,v,edata['weight']) for u,v,edata in G.edges(data=True) if 'weight' in edata ]
    [(1, 2, 4), (2, 3, 8)]

    **Reporting:**

    Simple graph information is obtained using methods.
    Iterator versions of many reporting methods exist for efficiency.
    Methods exist for reporting nodes(), edges(), neighbors() and degree()
    as well as the number of nodes and edges.

    For details on these and other miscellaneous methods, see below.
    """
    def __init__(self, data=None, **attr):
        """Initialize a graph with edges, name, graph attributes.

        Parameters
        ----------
        data : input graph
            Data to initialize graph.  If data=None (default) an empty
            graph is created.  The data can be an edge list, or any
            NetworkX graph object.  If the corresponding optional Python
            packages are installed the data can also be a NumPy matrix
            or 2d ndarray, a SciPy sparse matrix, or a PyGraphviz graph.
        name : string, optional (default='')
            An optional name for the graph.
        attr : keyword arguments, optional (default= no attributes)
            Attributes to add to graph as key=value pairs.

        See Also
        --------
        convert

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G = Graph(name='my graph')
        >>> e = [(1,2),(2,3),(3,4)] # list of edges
        >>> G = Graph(e)

        Arbitrary graph attribute pairs (key=value) may be assigned

        >>> G=Graph(e, day="Friday")
        >>> G.graph
        {'day': 'Friday'}

        """
        self.graph = {}   # dictionary for graph attributes
        self.node = {}    # empty node dict (created before convert)
        self.adj = {}     # empty adjacency dict
        # attempt to load graph with data
        if data is not None:
            to_networkx_graph(data,create_using=self)
        # load graph attributes (must be after convert)
        self.graph.update(attr)
        self.edge = self.adj

    @property
    def name(self):
        return self.graph.get('name','')
    
    @name.setter
    def name(self, s):
        self.graph['name']=s

    def __str__(self):
        """Return the graph name.

        Returns
        -------
        name : string
            The name of the graph.

        Examples
        --------
        >>> G = Graph(name='foo')
        >>> str(G)
        'foo'
        """
        return self.name

    def __iter__(self):
        """Iterate over the nodes. Use the expression 'for n in G'.

        Returns
        -------
        niter : iterator
            An iterator over all nodes in the graph.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        """
        return iter(self.adj)

    def __contains__(self,n):
        """Return True if n is a node, False otherwise. Use the expression
        'n in G'.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> 1 in G
        True
        """
        try:
            return n in self.adj
        except TypeError:
            return False

    def __len__(self):
        """Return the number of nodes. Use the expression 'len(G)'.

        Returns
        -------
        nnodes : int
            The number of nodes in the graph.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> len(G)
        4

        """
        return len(self.adj)

    def __getitem__(self, n):
        """Return a dict of neighbors of node n.  Use the expression 'G[n]'.

        Parameters
        ----------
        n : node
           A node in the graph.

        Returns
        -------
        adj_dict : dictionary
           The adjacency dictionary for nodes connected to n.

        Notes
        -----
        G[n] is similar to G.neighbors(n) but the internal data dictionary
        is returned instead of a list.

        Assigning G[n] will corrupt the internal graph data structure.
        Use G[n] for reading data only.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G[0]
        {1: {}}
        """
        return self.adj[n]


    def add_node(self, n, attr_dict=None, **attr):
        """Add a single node n and update node attributes.

        Parameters
        ----------
        n : node
            A node can be any hashable Python object except None.
        attr_dict : dictionary, optional (default= no attributes)
            Dictionary of node attributes.  Key/value pairs will
            update existing data associated with the node.
        attr : keyword arguments, optional
            Set or change attributes using key=value.

        See Also
        --------
        add_nodes_from

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_node(1)
        >>> G.add_node('Hello')
        >>> K3 = Graph([(0,1),(1,2),(2,0)])
        >>> G.add_node(K3)
        >>> G.number_of_nodes()
        3

        Use keywords set/change node attributes:

        >>> G.add_node(1,size=10)
        >>> G.add_node(3,weight=0.4,UTM=('13S',382871,3972649))

        Notes
        -----
        A hashable object is one that can be used as a key in a Python
        dictionary. This includes strings, numbers, tuples of strings
        and numbers, etc.

        On many platforms hashable items also include mutables such as
        NetworkX Graphs, though one should be careful that the hash
        doesn't change on mutables.
        """
        # set up attribute dict
        if attr_dict is None:
            attr_dict=attr
        else:
            try:
                attr_dict.update(attr)
            except AttributeError:
                raise NetworkXError(\
                    "The attr_dict argument must be a dictionary.")
        if n not in self.adj:
            self.adj[n] = {}
            self.node[n] = attr_dict
        else: # update attr even if node already exists
            self.node[n].update(attr_dict)


    def add_nodes_from(self, nodes, **attr):
        """Add multiple nodes.

        Parameters
        ----------
        nodes : iterable container
            A container of nodes (list, dict, set, etc.).
            OR
            A container of (node, attribute dict) tuples.
            Node attributes are updated using the attribute dict.
        attr : keyword arguments, optional (default= no attributes)
            Update attributes for all nodes in nodes.
            Node attributes specified in nodes as a tuple
            take precedence over attributes specified generally.

        See Also
        --------
        add_node

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_nodes_from('Hello')
        >>> K3 = Graph([(0,1),(1,2),(2,0)])
        >>> G.add_nodes_from(K3)
        >>> sorted(G.nodes(),key=str)
        [0, 1, 2, 'H', 'e', 'l', 'o']

        Use keywords to update specific node attributes for every node.

        >>> G.add_nodes_from([1,2], size=10)
        >>> G.add_nodes_from([3,4], weight=0.4)

        Use (node, attrdict) tuples to update attributes for specific
        nodes.

        >>> G.add_nodes_from([(1,dict(size=11)), (2,{'color':'blue'})])
        >>> G.node[1]['size']
        11
        >>> H = Graph()
        >>> H.add_nodes_from(G.nodes(data=True))
        >>> H.node[1]['size']
        11

        """
        for n in nodes:
            try:
                newnode=n not in self.adj
            except TypeError:
                nn,ndict = n
                if nn not in self.adj:
                    self.adj[nn] = {}
                    newdict = attr.copy()
                    newdict.update(ndict)
                    self.node[nn] = newdict
                else:
                    olddict = self.node[nn]
                    olddict.update(attr)
                    olddict.update(ndict)
                continue
            if newnode:
                self.adj[n] = {}
                self.node[n] = attr.copy()
            else:
                self.node[n].update(attr)

    def remove_node(self,n):
        """Remove node n.

        Removes the node n and all adjacent edges.
        Attempting to remove a non-existent node will raise an exception.

        Parameters
        ----------
        n : node
           A node in the graph

        Raises
        -------
        NetworkXError
           If n is not in the graph.

        See Also
        --------
        remove_nodes_from

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2])
        >>> G.edges()
        [(0, 1), (1, 2)]
        >>> G.remove_node(1)
        >>> G.edges()
        []

        """
        adj = self.adj
        try:
            nbrs = list(adj[n].keys()) # keys handles self-loops (allow mutation later)
            del self.node[n]
        except KeyError: # NetworkXError if n not in self
            raise NetworkXError("The node %s is not in the graph."%(n,))
        for u in nbrs:
            del adj[u][n]   # remove all edges n-u in graph
        del adj[n]          # now remove node


    def remove_nodes_from(self, nodes):
        """Remove multiple nodes.

        Parameters
        ----------
        nodes : iterable container
            A container of nodes (list, dict, set, etc.).  If a node
            in the container is not in the graph it is silently
            ignored.

        See Also
        --------
        remove_node

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2])
        >>> e = G.nodes()
        >>> e
        [0, 1, 2]
        >>> G.remove_nodes_from(e)
        >>> G.nodes()
        []

        """
        adj = self.adj
        for n in nodes:
            try:
                del self.node[n]
                for u in list(adj[n].keys()):   # keys() handles self-loops 
                    del adj[u][n]         #(allows mutation of dict in loop)
                del adj[n]
            except KeyError:
                pass


    def nodes_iter(self, data=False):
        """Return an iterator over the nodes.

        Parameters
        ----------
        data : boolean, optional (default=False)
               If False the iterator returns nodes.  If True
               return a two-tuple of node and node data dictionary

        Returns
        -------
        niter : iterator
            An iterator over nodes.  If data=True the iterator gives
            two-tuples containing (node, node data, dictionary)

        Notes
        -----
        If the node data is not required it is simpler and equivalent
        to use the expression 'for n in G'.

        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2])

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2])

        >>> [d for n,d in G.nodes_iter(data=True)]
        [{}, {}, {}]
        """
        if data:
            return iter(self.node.items())
        return iter(self.adj)

    def nodes(self, data=False):
        """Return a list of the nodes in the graph.

        Parameters
        ----------
        data : boolean, optional (default=False)
               If False return a list of nodes.  If True return a
               two-tuple of node and node data dictionary

        Returns
        -------
        nlist : list
            A list of nodes.  If data=True a list of two-tuples containing
            (node, node data dictionary).

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2])
        >>> G.nodes()
        [0, 1, 2]
        >>> G.add_node(1, time='5pm')
        >>> G.nodes(data=True)
        [(0, {}), (1, {'time': '5pm'}), (2, {})]
        """
        return list(self.nodes_iter(data=data))

    def number_of_nodes(self):
        """Return the number of nodes in the graph.

        Returns
        -------
        nnodes : int
            The number of nodes in the graph.

        See Also
        --------
        order, __len__  which are identical

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2])
        >>> len(G)
        3
        """
        return len(self.adj)

    def order(self):
        """Return the number of nodes in the graph.

        Returns
        -------
        nnodes : int
            The number of nodes in the graph.

        See Also
        --------
        number_of_nodes, __len__  which are identical

        """
        return len(self.adj)

    def has_node(self, n):
        """Return True if the graph contains the node n.

        Parameters
        ----------
        n : node

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2])
        >>> G.has_node(0)
        True

        It is more readable and simpler to use

        >>> 0 in G
        True

        """
        try:
            return n in self.adj
        except TypeError:
            return False

    def add_edge(self, u, v, attr_dict=None, **attr):
        """Add an edge between u and v.

        The nodes u and v will be automatically added if they are
        not already in the graph.

        Edge attributes can be specified with keywords or by providing
        a dictionary with key/value pairs.  See examples below.

        Parameters
        ----------
        u,v : nodes
            Nodes can be, for example, strings or numbers.
            Nodes must be hashable (and not None) Python objects.
        attr_dict : dictionary, optional (default= no attributes)
            Dictionary of edge attributes.  Key/value pairs will
            update existing data associated with the edge.
        attr : keyword arguments, optional
            Edge data (or labels or objects) can be assigned using
            keyword arguments.

        See Also
        --------
        add_edges_from : add a collection of edges

        Notes
        -----
        Adding an edge that already exists updates the edge data.

        Many NetworkX algorithms designed for weighted graphs use as
        the edge weight a numerical value assigned to a keyword
        which by default is 'weight'.

        Examples
        --------
        The following all add the edge e=(1,2) to graph G:

        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> e = (1,2)
        >>> G.add_edge(1, 2)           # explicit two-node form
        >>> G.add_edge(*e)             # single edge as tuple of two nodes
        >>> G.add_edges_from( [(1,2)] ) # add edges from iterable container

        Associate data to edges using keywords:

        >>> G.add_edge(1, 2, weight=3)
        >>> G.add_edge(1, 3, weight=7, capacity=15, length=342.7)
        """
        # set up attribute dictionary
        if attr_dict is None:
            attr_dict=attr
        else:
            try:
                attr_dict.update(attr)
            except AttributeError:
                raise NetworkXError(\
                    "The attr_dict argument must be a dictionary.")
        # add nodes
        if u not in self.adj:
            self.adj[u] = {}
            self.node[u] = {}
        if v not in self.adj:
            self.adj[v] = {}
            self.node[v] = {}
        # add the edge
        datadict=self.adj[u].get(v,{})
        datadict.update(attr_dict)
        self.adj[u][v] = datadict
        self.adj[v][u] = datadict


    def add_edges_from(self, ebunch, attr_dict=None, **attr):
        """Add all the edges in ebunch.

        Parameters
        ----------
        ebunch : container of edges
            Each edge given in the container will be added to the
            graph. The edges must be given as as 2-tuples (u,v) or
            3-tuples (u,v,d) where d is a dictionary containing edge
            data.
        attr_dict : dictionary, optional (default= no attributes)
            Dictionary of edge attributes.  Key/value pairs will
            update existing data associated with each edge.
        attr : keyword arguments, optional
            Edge data (or labels or objects) can be assigned using
            keyword arguments.


        See Also
        --------
        add_edge : add a single edge
        add_weighted_edges_from : convenient way to add weighted edges

        Notes
        -----
        Adding the same edge twice has no effect but any edge data
        will be updated when each duplicate edge is added.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edges_from([(0,1),(1,2)]) # using a list of edge tuples
        >>> e = zip(range(0,3),range(1,4))
        >>> G.add_edges_from(e) # Add the path graph 0-1-2-3

        Associate data to edges

        >>> G.add_edges_from([(1,2),(2,3)], weight=3)
        >>> G.add_edges_from([(3,4),(1,4)], label='WN2898')
        """
        # set up attribute dict
        if attr_dict is None:
            attr_dict=attr
        else:
            try:
                attr_dict.update(attr)
            except AttributeError:
                raise NetworkXError(\
                    "The attr_dict argument must be a dictionary.")
        # process ebunch
        for e in ebunch:
            ne=len(e)
            if ne==3:
                u,v,dd = e
            elif ne==2:
                u,v = e
                dd = {}
            else:
                raise NetworkXError(\
                    "Edge tuple %s must be a 2-tuple or 3-tuple."%(e,))
            if u not in self.adj:
                self.adj[u] = {}
                self.node[u] = {}
            if v not in self.adj:
                self.adj[v] = {}
                self.node[v] = {}
            datadict=self.adj[u].get(v,{})
            datadict.update(attr_dict)
            datadict.update(dd)
            self.adj[u][v] = datadict
            self.adj[v][u] = datadict


    def add_weighted_edges_from(self, ebunch, weight='weight', **attr):
        """Add all the edges in ebunch as weighted edges with specified
        weights.

        Parameters
        ----------
        ebunch : container of edges
            Each edge given in the list or container will be added
            to the graph. The edges must be given as 3-tuples (u,v,w)
            where w is a number.
        weight : string, optional (default= 'weight')
            The attribute name for the edge weights to be added.
        attr : keyword arguments, optional (default= no attributes)
            Edge attributes to add/update for all edges.

        See Also
        --------
        add_edge : add a single edge
        add_edges_from : add multiple edges

        Notes
        -----
        Adding the same edge twice for Graph/DiGraph simply updates 
        the edge data.  For MultiGraph/MultiDiGraph, duplicate edges 
        are stored.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_weighted_edges_from([(0,1,3.0),(1,2,7.5)])
        """
        self.add_edges_from(((u,v,{weight:d}) for u,v,d in ebunch),**attr)

    def remove_edge(self, u, v):
        """Remove the edge between u and v.

        Parameters
        ----------
        u,v: nodes
            Remove the edge between nodes u and v.

        Raises
        ------
        NetworkXError
            If there is not an edge between u and v.

        See Also
        --------
        remove_edges_from : remove a collection of edges

        Examples
        --------
        >>> G = Graph()   # or DiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.remove_edge(0,1)
        >>> e = (1,2)
        >>> G.remove_edge(*e) # unpacks e from an edge tuple
        >>> e = (2,3,{'weight':7}) # an edge with attribute data
        >>> G.remove_edge(*e[:2]) # select first part of edge tuple
        """
        try:
            del self.adj[u][v]
            if u != v:  # self-loop needs only one entry removed
                del self.adj[v][u]
        except KeyError:
            raise NetworkXError("The edge %s-%s is not in the graph"%(u,v))



    def remove_edges_from(self, ebunch):
        """Remove all edges specified in ebunch.

        Parameters
        ----------
        ebunch: list or container of edge tuples
            Each edge given in the list or container will be removed
            from the graph. The edges can be:

                - 2-tuples (u,v) edge between u and v.
                - 3-tuples (u,v,k) where k is ignored.

        See Also
        --------
        remove_edge : remove a single edge

        Notes
        -----
        Will fail silently if an edge in ebunch is not in the graph.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> ebunch=[(1,2),(2,3)]
        >>> G.remove_edges_from(ebunch)
        """
        for e in ebunch:
            u,v = e[:2]  # ignore edge data if present
            if u in self.adj and v in self.adj[u]:
                del self.adj[u][v]
                if u != v:  # self loop needs only one entry removed
                    del self.adj[v][u]


    def has_edge(self, u, v):
        """Return True if the edge (u,v) is in the graph.

        Parameters
        ----------
        u,v : nodes
            Nodes can be, for example, strings or numbers.
            Nodes must be hashable (and not None) Python objects.

        Returns
        -------
        edge_ind : bool
            True if edge is in the graph, False otherwise.

        Examples
        --------
        Can be called either using two nodes u,v or edge tuple (u,v)

        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.has_edge(0,1)  # using two nodes
        True
        >>> e = (0,1)
        >>> G.has_edge(*e)  #  e is a 2-tuple (u,v)
        True
        >>> e = (0,1,{'weight':7})
        >>> G.has_edge(*e[:2])  # e is a 3-tuple (u,v,data_dictionary)
        True

        The following syntax are all equivalent:

        >>> G.has_edge(0,1)
        True
        >>> 1 in G[0]  # though this gives KeyError if 0 not in G
        True

        """
        try:
            return v in self.adj[u]
        except KeyError:
            return False


    def neighbors(self, n):
        """Return a list of the nodes connected to the node n.

        Parameters
        ----------
        n : node
           A node in the graph

        Returns
        -------
        nlist : list
            A list of nodes that are adjacent to n.

        Raises
        ------
        NetworkXError
            If the node n is not in the graph.

        Notes
        -----
        It is usually more convenient (and faster) to access the
        adjacency dictionary as G[n]:

        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edge('a','b',weight=7)
        >>> G['a']
        {'b': {'weight': 7}}

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.neighbors(0)
        [1]

        """
        try:
            return list(self.adj[n])
        except KeyError:
            raise NetworkXError("The node %s is not in the graph."%(n,))

    def neighbors_iter(self, n):
        """Return an iterator over all neighbors of node n.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> [n for n in G.neighbors_iter(0)]
        [1]

        Notes
        -----
        It is faster to use the idiom "in G[0]", e.g.

        >>> G = path_graph(4)
        >>> [n for n in G[0]]
        [1]
        """
        try:
            return iter(self.adj[n])
        except KeyError:
            raise NetworkXError("The node %s is not in the graph."%(n,))

    def edges(self, nbunch=None, data=False):
        """Return a list of edges.

        Edges are returned as tuples with optional data
        in the order (node, neighbor, data).

        Parameters
        ----------
        nbunch : iterable container, optional (default= all nodes)
            A container of nodes.  The container will be iterated
            through once.
        data : bool, optional (default=False)
            Return two tuples (u,v) (False) or three-tuples (u,v,data) (True).

        Returns
        --------
        edge_list: list of edge tuples
            Edges that are adjacent to any node in nbunch, or a list
            of all edges if nbunch is not specified.

        See Also
        --------
        edges_iter : return an iterator over the edges

        Notes
        -----
        Nodes in nbunch that are not in the graph will be (quietly) ignored.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.edges()
        [(0, 1), (1, 2), (2, 3)]
        >>> G.edges(data=True) # default edge data is {} (empty dictionary)
        [(0, 1, {}), (1, 2, {}), (2, 3, {})]
        >>> G.edges([0,3])
        [(0, 1), (3, 2)]
        >>> G.edges(0)
        [(0, 1)]

        """
        return list(self.edges_iter(nbunch, data))

    def edges_iter(self, nbunch=None, data=False):
        """Return an iterator over the edges.

        Edges are returned as tuples with optional data
        in the order (node, neighbor, data).

        Parameters
        ----------
        nbunch : iterable container, optional (default= all nodes)
            A container of nodes.  The container will be iterated
            through once.
        data : bool, optional (default=False)
            If True, return edge attribute dict in 3-tuple (u,v,data).

        Returns
        -------
        edge_iter : iterator
            An iterator of (u,v) or (u,v,d) tuples of edges.

        See Also
        --------
        edges : return a list of edges

        Notes
        -----
        Nodes in nbunch that are not in the graph will be (quietly) ignored.

        Examples
        --------
        >>> G = Graph()   # or MultiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> [e for e in G.edges_iter()]
        [(0, 1), (1, 2), (2, 3)]
        >>> list(G.edges_iter(data=True)) # default data is {} (empty dict)
        [(0, 1, {}), (1, 2, {}), (2, 3, {})]
        >>> list(G.edges_iter([0,3]))
        [(0, 1), (3, 2)]
        >>> list(G.edges_iter(0))
        [(0, 1)]

        """
        seen={}     # helper dict to keep track of multiply stored edges
        if nbunch is None:
            nodes_nbrs = self.adj.items()
        else:
            nodes_nbrs=((n,self.adj[n]) for n in self.nbunch_iter(nbunch))
        if data:
            for n,nbrs in nodes_nbrs:
                for nbr,data in nbrs.items():
                    if nbr not in seen:
                        yield (n,nbr,data)
                seen[n]=1
        else:
            for n,nbrs in nodes_nbrs:
                for nbr in nbrs:
                    if nbr not in seen:
                        yield (n,nbr)
                seen[n] = 1
        del seen


    def get_edge_data(self, u, v, default=None):
        """Return the attribute dictionary associated with edge (u,v).

        Parameters
        ----------
        u,v : nodes
        default:  any Python object (default=None)
            Value to return if the edge (u,v) is not found.

        Returns
        -------
        edge_dict : dictionary
            The edge attribute dictionary.

        Notes
        -----
        It is faster to use G[u][v].

        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G[0][1]
        {}

        Warning: Assigning G[u][v] corrupts the graph data structure.
        But it is safe to assign attributes to that dictionary,

        >>> G[0][1]['weight'] = 7
        >>> G[0][1]['weight']
        7
        >>> G[1][0]['weight']
        7

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.get_edge_data(0,1) # default edge data is {}
        {}
        >>> e = (0,1)
        >>> G.get_edge_data(*e) # tuple form
        {}
        >>> G.get_edge_data('a','b',default=0) # edge not in graph, return 0
        0
        """
        try:
            return self.adj[u][v]
        except KeyError:
            return default

    def adjacency_list(self):
        """Return an adjacency list representation of the graph.

        The output adjacency list is in the order of G.nodes().
        For directed graphs, only outgoing adjacencies are included.

        Returns
        -------
        adj_list : lists of lists
            The adjacency structure of the graph as a list of lists.

        See Also
        --------
        adjacency_iter

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.adjacency_list() # in order given by G.nodes()
        [[1], [0, 2], [1, 3], [2]]

        """
        return list(map(list,iter(self.adj.values())))

    def adjacency_iter(self):
        """Return an iterator of (node, adjacency dict) tuples for all nodes.

        This is the fastest way to look at every edge.
        For directed graphs, only outgoing adjacencies are included.

        Returns
        -------
        adj_iter : iterator
           An iterator of (node, adjacency dictionary) for all nodes in
           the graph.

        See Also
        --------
        adjacency_list

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> [(n,nbrdict) for n,nbrdict in G.adjacency_iter()]
        [(0, {1: {}}), (1, {0: {}, 2: {}}), (2, {1: {}, 3: {}}), (3, {2: {}})]

        """
        return iter(self.adj.items())

    def degree(self, nbunch=None, weight=None):
        """Return the degree of a node or nodes.

        The node degree is the number of edges adjacent to that node.

        Parameters
        ----------
        nbunch : iterable container, optional (default=all nodes)
            A container of nodes.  The container will be iterated
            through once.

        weight : string or None, optional (default=None)
           The edge attribute that holds the numerical value used 
           as a weight.  If None, then each edge has weight 1.
           The degree is the sum of the edge weights adjacent to the node.

        Returns
        -------
        nd : dictionary, or number
            A dictionary with nodes as keys and degree as values or
            a number if a single node is specified.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.degree(0)
        1
        >>> G.degree([0,1])
        {0: 1, 1: 2}
        >>> list(G.degree([0,1]).values())
        [1, 2]

        """
        if nbunch in self:      # return a single node
            return next(self.degree_iter(nbunch,weight))[1]
        else:           # return a dict
            return dict(self.degree_iter(nbunch,weight))

    def degree_iter(self, nbunch=None, weight=None):
        """Return an iterator for (node, degree).

        The node degree is the number of edges adjacent to the node.

        Parameters
        ----------
        nbunch : iterable container, optional (default=all nodes)
            A container of nodes.  The container will be iterated
            through once.

        weight : string or None, optional (default=None)
           The edge attribute that holds the numerical value used 
           as a weight.  If None, then each edge has weight 1.
           The degree is the sum of the edge weights adjacent to the node.

        Returns
        -------
        nd_iter : an iterator
            The iterator returns two-tuples of (node, degree).

        See Also
        --------
        degree

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> list(G.degree_iter(0)) # node 0 with degree 1
        [(0, 1)]
        >>> list(G.degree_iter([0,1]))
        [(0, 1), (1, 2)]

        """
        if nbunch is None:
            nodes_nbrs = self.adj.items()
        else:
            nodes_nbrs=((n,self.adj[n]) for n in self.nbunch_iter(nbunch))
  
        if weight is None:
            for n,nbrs in nodes_nbrs:
                yield (n,len(nbrs)+(n in nbrs)) # return tuple (n,degree)
        else:
        # edge weighted graph - degree is sum of nbr edge weights
            for n,nbrs in nodes_nbrs:
                yield (n, sum((nbrs[nbr].get(weight,1) for nbr in nbrs)) +
                              (n in nbrs and nbrs[n].get(weight,1)))


    def clear(self):
        """Remove all nodes and edges from the graph.

        This also removes the name, and all graph, node, and edge attributes.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.clear()
        >>> G.nodes()
        []
        >>> G.edges()
        []

        """
        self.name = ''
        self.adj.clear()
        self.node.clear()
        self.graph.clear()

    def copy(self):
        """Return a copy of the graph.

        Returns
        -------
        G : Graph
            A copy of the graph.

        See Also
        --------
        to_directed: return a directed copy of the graph.

        Notes
        -----
        This makes a complete copy of the graph including all of the
        node or edge attributes.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> H = G.copy()

        """
        return deepcopy(self)

    def is_multigraph(self):
        """Return True if graph is a multigraph, False otherwise."""
        return False


    def is_directed(self):
        """Return True if graph is directed, False otherwise."""
        return False

    def to_directed(self):
        """Return a directed representation of the graph.

        Returns
        -------
        G : DiGraph
            A directed graph with the same name, same nodes, and with
            each edge (u,v,data) replaced by two directed edges
            (u,v,data) and (v,u,data).

        Notes
        -----
        This returns a "deepcopy" of the edge, node, and
        graph attributes which attempts to completely copy
        all of the data and references.

        This is in contrast to the similar D=DiGraph(G) which returns a
        shallow copy of the data.

        See the Python copy module for more information on shallow
        and deep copies, http://docs.python.org/library/copy.html.

        Examples
        --------
        >>> G = Graph()   # or MultiGraph, etc
        >>> G.add_path([0,1])
        >>> H = G.to_directed()
        >>> H.edges()
        [(0, 1), (1, 0)]

        If already directed, return a (deep) copy

        >>> G = DiGraph()   # or MultiDiGraph, etc
        >>> G.add_path([0,1])
        >>> H = G.to_directed()
        >>> H.edges()
        [(0, 1)]
        """
        G=DiGraph()
        G.name=self.name
        G.add_nodes_from(self)
        G.add_edges_from( ((u,v,deepcopy(data)) 
                           for u,nbrs in self.adjacency_iter() 
                           for v,data in nbrs.items()) )
        G.graph=deepcopy(self.graph)
        G.node=deepcopy(self.node)
        return G

    def to_undirected(self):
        """Return an undirected copy of the graph.

        Returns
        -------
        G : Graph/MultiGraph
            A deepcopy of the graph.

        See Also
        --------
        copy, add_edge, add_edges_from

        Notes
        -----
        This returns a "deepcopy" of the edge, node, and
        graph attributes which attempts to completely copy
        all of the data and references.

        This is in contrast to the similar G=DiGraph(D) which returns a
        shallow copy of the data.

        See the Python copy module for more information on shallow
        and deep copies, http://docs.python.org/library/copy.html.

        Examples
        --------
        >>> G = Graph()   # or MultiGraph, etc
        >>> G.add_path([0,1])
        >>> H = G.to_directed()
        >>> H.edges()
        [(0, 1), (1, 0)]
        >>> G2 = H.to_undirected()
        >>> G2.edges()
        [(0, 1)]
        """
        return deepcopy(self)

    def subgraph(self, nbunch):
        """Return the subgraph induced on nodes in nbunch.

        The induced subgraph of the graph contains the nodes in nbunch
        and the edges between those nodes.

        Parameters
        ----------
        nbunch : list, iterable
            A container of nodes which will be iterated through once.

        Returns
        -------
        G : Graph
            A subgraph of the graph with the same edge attributes.

        Notes
        -----
        The graph, edge or node attributes just point to the original graph.
        So changes to the node or edge structure will not be reflected in
        the original graph while changes to the attributes will.

        To create a subgraph with its own copy of the edge/node attributes use:
        Graph(G.subgraph(nbunch))

        If edge attributes are containers, a deep copy can be obtained using:
        G.subgraph(nbunch).copy()

        For an inplace reduction of a graph to a subgraph you can remove nodes:
        G.remove_nodes_from([ n in G if n not in set(nbunch)])

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> H = G.subgraph([0,1,2])
        >>> H.edges()
        [(0, 1), (1, 2)]
        """
        bunch =self.nbunch_iter(nbunch)
        # create new graph and copy subgraph into it
        H = self.__class__()
        # namespace shortcuts for speed
        H_adj=H.adj
        self_adj=self.adj
        # add nodes and edges (undirected method)
        for n in bunch:
            Hnbrs={}
            H_adj[n]=Hnbrs
            for nbr,d in self_adj[n].items():
                if nbr in H_adj:
                    # add both representations of edge: n-nbr and nbr-n
                    Hnbrs[nbr]=d
                    H_adj[nbr][n]=d
        # copy node and attribute dictionaries
        for n in H:
            H.node[n]=self.node[n]
        H.graph=self.graph
        return H


    def nodes_with_selfloops(self):
        """Return a list of nodes with self loops.

        A node with a self loop has an edge with both ends adjacent
        to that node.

        Returns
        -------
        nodelist : list
            A list of nodes with self loops.

        See Also
        --------
        selfloop_edges, number_of_selfloops

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edge(1,1)
        >>> G.add_edge(1,2)
        >>> G.nodes_with_selfloops()
        [1]
        """
        return [ n for n,nbrs in self.adj.items() if n in nbrs ]

    def selfloop_edges(self, data=False):
        """Return a list of selfloop edges.

        A selfloop edge has the same node at both ends.

        Parameters
        -----------
        data : bool, optional (default=False)
            Return selfloop edges as two tuples (u,v) (data=False)
            or three-tuples (u,v,data) (data=True)

        Returns
        -------
        edgelist : list of edge tuples
            A list of all selfloop edges.

        See Also
        --------
        nodes_with_selfloops, number_of_selfloops

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edge(1,1)
        >>> G.add_edge(1,2)
        >>> G.selfloop_edges()
        [(1, 1)]
        >>> G.selfloop_edges(data=True)
        [(1, 1, {})]
        """
        if data:
            return [ (n,n,nbrs[n])
                     for n,nbrs in self.adj.items() if n in nbrs ]
        else:
            return [ (n,n)
                     for n,nbrs in self.adj.items() if n in nbrs ]


    def number_of_selfloops(self):
        """Return the number of selfloop edges.

        A selfloop edge has the same node at both ends.

        Returns
        -------
        nloops : int
            The number of selfloops.

        See Also
        --------
        nodes_with_selfloops, selfloop_edges

        Examples
        --------
        >>> G=Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edge(1,1)
        >>> G.add_edge(1,2)
        >>> G.number_of_selfloops()
        1
        """
        return len(self.selfloop_edges())


    def size(self, weight=None):
        """Return the number of edges.

        Parameters
        ----------
        weight : string or None, optional (default=None)
           The edge attribute that holds the numerical value used 
           as a weight.  If None, then each edge has weight 1.

        Returns
        -------
        nedges : int
            The number of edges of sum of edge weights in the graph.

        See Also
        --------
        number_of_edges

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.size()
        3

        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edge('a','b',weight=2)
        >>> G.add_edge('b','c',weight=4)
        >>> G.size()
        2
        >>> G.size(weight='weight')
        6.0
        """
        s=sum(self.degree(weight=weight).values())/2
        if weight is None:
            return int(s)
        else:
            return float(s)

    def number_of_edges(self, u=None, v=None):
        """Return the number of edges between two nodes.

        Parameters
        ----------
        u,v : nodes, optional (default=all edges)
            If u and v are specified, return the number of edges between
            u and v. Otherwise return the total number of all edges.

        Returns
        -------
        nedges : int
            The number of edges in the graph.  If nodes u and v are specified
            return the number of edges between those nodes.

        See Also
        --------
        size

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.number_of_edges()
        3
        >>> G.number_of_edges(0,1)
        1
        >>> e = (0,1)
        >>> G.number_of_edges(*e)
        1
        """
        if u is None: return int(self.size())
        if v in self.adj[u]:
            return 1
        else:
            return 0


    def add_star(self, nodes, **attr):
        """Add a star.

        The first node in nodes is the middle of the star.  It is connected
        to all other nodes.

        Parameters
        ----------
        nodes : iterable container
            A container of nodes.
        attr : keyword arguments, optional (default= no attributes)
            Attributes to add to every edge in star.

        See Also
        --------
        add_path, add_cycle

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_star([0,1,2,3])
        >>> G.add_star([10,11,12],weight=2)

        """
        nlist = list(nodes)
        v=nlist[0]
        edges=((v,n) for n in nlist[1:])
        self.add_edges_from(edges, **attr)

    def add_path(self, nodes, **attr):
        """Add a path.

        Parameters
        ----------
        nodes : iterable container
            A container of nodes.  A path will be constructed from
            the nodes (in order) and added to the graph.
        attr : keyword arguments, optional (default= no attributes)
            Attributes to add to every edge in path.

        See Also
        --------
        add_star, add_cycle

        Examples
        --------
        >>> G=Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.add_path([10,11,12],weight=7)

        """
        nlist = list(nodes)
        edges=zip(nlist[:-1],nlist[1:])
        self.add_edges_from(edges, **attr)

    def add_cycle(self, nodes, **attr):
        """Add a cycle.

        Parameters
        ----------
        nodes: iterable container
            A container of nodes.  A cycle will be constructed from
            the nodes (in order) and added to the graph.
        attr : keyword arguments, optional (default= no attributes)
            Attributes to add to every edge in cycle.

        See Also
        --------
        add_path, add_star

        Examples
        --------
        >>> G=Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_cycle([0,1,2,3])
        >>> G.add_cycle([10,11,12],weight=7)

        """
        nlist = list(nodes)
        edges=zip(nlist,nlist[1:]+[nlist[0]])
        self.add_edges_from(edges, **attr)


    def nbunch_iter(self, nbunch=None):
        """Return an iterator of nodes contained in nbunch that are
        also in the graph.

        The nodes in nbunch are checked for membership in the graph
        and if not are silently ignored.

        Parameters
        ----------
        nbunch : iterable container, optional (default=all nodes)
            A container of nodes.  The container will be iterated
            through once.

        Returns
        -------
        niter : iterator
            An iterator over nodes in nbunch that are also in the graph.
            If nbunch is None, iterate over all nodes in the graph.

        Raises
        ------
        NetworkXError
            If nbunch is not a node or or sequence of nodes.
            If a node in nbunch is not hashable.

        See Also
        --------
        Graph.__iter__

        Notes
        -----
        When nbunch is an iterator, the returned iterator yields values
        directly from nbunch, becoming exhausted when nbunch is exhausted.

        To test whether nbunch is a single node, one can use
        "if nbunch in self:", even after processing with this routine.

        If nbunch is not a node or a (possibly empty) sequence/iterator
        or None, a NetworkXError is raised.  Also, if any object in
        nbunch is not hashable, a NetworkXError is raised.
        """
        if nbunch is None:   # include all nodes via iterator
            bunch=iter(self.adj.keys())
        elif nbunch in self: # if nbunch is a single node
            bunch=iter([nbunch])
        else:                # if nbunch is a sequence of nodes
            def bunch_iter(nlist,adj):
                try:
                    for n in nlist:
                        if n in adj:
                            yield n
                except TypeError as e:
                    message=e.args[0]
#                    sys.stdout.write(message)
                    # capture error for non-sequence/iterator nbunch.
                    if 'iter' in message:
                        raise NetworkXError(\
                            "nbunch is not a node or a sequence of nodes.")
                    # capture error for unhashable node.
                    elif 'hashable' in message:
                        raise NetworkXError(\
                            "Node %s in the sequence nbunch is not a valid node."%n)
                    else: 
                        raise 
            bunch=bunch_iter(nbunch,self.adj)
        return bunch

'''
============================================================
DiGraph, from networkx
============================================================
'''
"""Base class for directed graphs."""
#    Copyright (C) 2004-2011 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

class DiGraph(Graph):
    """
    Base class for directed graphs.

    A DiGraph stores nodes and edges with optional data, or attributes.

    DiGraphs hold directed edges.  Self loops are allowed but multiple
    (parallel) edges are not.

    Nodes can be arbitrary (hashable) Python objects with optional
    key/value attributes.

    Edges are represented as links between nodes with optional
    key/value attributes.

    Parameters
    ----------
    data : input graph
        Data to initialize graph.  If data=None (default) an empty
        graph is created.  The data can be an edge list, or any
        NetworkX graph object.  If the corresponding optional Python
        packages are installed the data can also be a NumPy matrix
        or 2d ndarray, a SciPy sparse matrix, or a PyGraphviz graph.
    attr : keyword arguments, optional (default= no attributes)
        Attributes to add to graph as key=value pairs.

    See Also
    --------
    Graph
    MultiGraph
    MultiDiGraph

    Examples
    --------
    Create an empty graph structure (a "null graph") with no nodes and
    no edges.

    >>> G = DiGraph()

    G can be grown in several ways.

    **Nodes:**

    Add one node at a time:

    >>> G.add_node(1)

    Add the nodes from any container (a list, dict, set or
    even the lines from a file or the nodes from another graph).

    >>> G.add_nodes_from([2,3])
    >>> G.add_nodes_from(range(100,110))
    >>> H=Graph()
    >>> H.add_path([0,1,2,3,4,5,6,7,8,9])
    >>> G.add_nodes_from(H)

    In addition to strings and integers any hashable Python object
    (except None) can represent a node, e.g. a customized node object,
    or even another Graph.

    >>> G.add_node(H)

    **Edges:**

    G can also be grown by adding edges.

    Add one edge,

    >>> G.add_edge(1, 2)

    a list of edges,

    >>> G.add_edges_from([(1,2),(1,3)])

    or a collection of edges,

    >>> G.add_edges_from(H.edges())

    If some edges connect nodes not yet in the graph, the nodes
    are added automatically.  There are no errors when adding
    nodes or edges that already exist.

    **Attributes:**

    Each graph, node, and edge can hold key/value attribute pairs
    in an associated attribute dictionary (the keys must be hashable).
    By default these are empty, but can be added or changed using
    add_edge, add_node or direct manipulation of the attribute
    dictionaries named graph, node and edge respectively.

    >>> G = DiGraph(day="Friday")
    >>> G.graph
    {'day': 'Friday'}

    Add node attributes using add_node(), add_nodes_from() or G.node

    >>> G.add_node(1, time='5pm')
    >>> G.add_nodes_from([3], time='2pm')
    >>> G.node[1]
    {'time': '5pm'}
    >>> G.node[1]['room'] = 714
    >>> G.nodes(data=True)
    [(1, {'room': 714, 'time': '5pm'}), (3, {'time': '2pm'})]

    Warning: adding a node to G.node does not add it to the graph.

    Add edge attributes using add_edge(), add_edges_from(), subscript
    notation, or G.edge.

    >>> G.add_edge(1, 2, weight=4.7 )
    >>> G.add_edges_from([(3,4),(4,5)], color='red')
    >>> G.add_edges_from([(1,2,{'color':'blue'}), (2,3,{'weight':8})])
    >>> G[1][2]['weight'] = 4.7
    >>> G.edge[1][2]['weight'] = 4

    **Shortcuts:**

    Many common graph features allow python syntax to speed reporting.

    >>> 1 in G     # check if node in graph
    True
    >>> [n for n in G if n<3]   # iterate through nodes
    [1, 2]
    >>> len(G)  # number of nodes in graph
    5
    >>> G[1] # adjacency dict keyed by neighbor to edge attributes
    ...            # Note: you should not change this dict manually!
    {2: {'color': 'blue', 'weight': 4}}

    The fastest way to traverse all edges of a graph is via
    adjacency_iter(), but the edges() method is often more convenient.

    >>> for n,nbrsdict in G.adjacency_iter():
    ...     for nbr,eattr in nbrsdict.items():
    ...        if 'weight' in eattr:
    ...            (n,nbr,eattr['weight'])
    (1, 2, 4)
    (2, 3, 8)
    >>> [ (u,v,edata['weight']) for u,v,edata in G.edges(data=True) if 'weight' in edata ]
    [(1, 2, 4), (2, 3, 8)]

    **Reporting:**

    Simple graph information is obtained using methods.
    Iterator versions of many reporting methods exist for efficiency.
    Methods exist for reporting nodes(), edges(), neighbors() and degree()
    as well as the number of nodes and edges.

    For details on these and other miscellaneous methods, see below.
    """
    def __init__(self, data=None, **attr):
        """Initialize a graph with edges, name, graph attributes.

        Parameters
        ----------
        data : input graph
            Data to initialize graph.  If data=None (default) an empty
            graph is created.  The data can be an edge list, or any
            NetworkX graph object.  If the corresponding optional Python
            packages are installed the data can also be a NumPy matrix
            or 2d ndarray, a SciPy sparse matrix, or a PyGraphviz graph.
        name : string, optional (default='')
            An optional name for the graph.
        attr : keyword arguments, optional (default= no attributes)
            Attributes to add to graph as key=value pairs.

        See Also
        --------
        convert

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G = Graph(name='my graph')
        >>> e = [(1,2),(2,3),(3,4)] # list of edges
        >>> G = Graph(e)

        Arbitrary graph attribute pairs (key=value) may be assigned

        >>> G=Graph(e, day="Friday")
        >>> G.graph
        {'day': 'Friday'}

        """
        self.graph = {} # dictionary for graph attributes
        self.node = {} # dictionary for node attributes
        # We store two adjacency lists:
        # the  predecessors of node n are stored in the dict self.pred
        # the successors of node n are stored in the dict self.succ=self.adj
        self.adj = {}  # empty adjacency dictionary
        self.pred = {}  # predecessor
        self.succ = self.adj  # successor

        # attempt to load graph with data
        if data is not None:
            to_networkx_graph(data,create_using=self)
        # load graph attributes (must be after convert)
        self.graph.update(attr)
        self.edge=self.adj


    def add_node(self, n, attr_dict=None, **attr):
        """Add a single node n and update node attributes.

        Parameters
        ----------
        n : node
            A node can be any hashable Python object except None.
        attr_dict : dictionary, optional (default= no attributes)
            Dictionary of node attributes.  Key/value pairs will
            update existing data associated with the node.
        attr : keyword arguments, optional
            Set or change attributes using key=value.

        See Also
        --------
        add_nodes_from

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_node(1)
        >>> G.add_node('Hello')
        >>> K3 = Graph([(0,1),(1,2),(2,0)])
        >>> G.add_node(K3)
        >>> G.number_of_nodes()
        3

        Use keywords set/change node attributes:

        >>> G.add_node(1,size=10)
        >>> G.add_node(3,weight=0.4,UTM=('13S',382871,3972649))

        Notes
        -----
        A hashable object is one that can be used as a key in a Python
        dictionary. This includes strings, numbers, tuples of strings
        and numbers, etc.

        On many platforms hashable items also include mutables such as
        NetworkX Graphs, though one should be careful that the hash
        doesn't change on mutables.
        """
        # set up attribute dict
        if attr_dict is None:
            attr_dict=attr
        else:
            try:
                attr_dict.update(attr)
            except AttributeError:
                raise NetworkXError(\
                    "The attr_dict argument must be a dictionary.")
        if n not in self.succ:
            self.succ[n] = {}
            self.pred[n] = {}
            self.node[n] = attr_dict
        else: # update attr even if node already exists
            self.node[n].update(attr_dict)


    def add_nodes_from(self, nodes, **attr):
        """Add multiple nodes.

        Parameters
        ----------
        nodes : iterable container
            A container of nodes (list, dict, set, etc.).
            OR
            A container of (node, attribute dict) tuples.
            Node attributes are updated using the attribute dict.
        attr : keyword arguments, optional (default= no attributes)
            Update attributes for all nodes in nodes.
            Node attributes specified in nodes as a tuple
            take precedence over attributes specified generally.

        See Also
        --------
        add_node

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_nodes_from('Hello')
        >>> K3 = Graph([(0,1),(1,2),(2,0)])
        >>> G.add_nodes_from(K3)
        >>> sorted(G.nodes(),key=str)
        [0, 1, 2, 'H', 'e', 'l', 'o']

        Use keywords to update specific node attributes for every node.

        >>> G.add_nodes_from([1,2], size=10)
        >>> G.add_nodes_from([3,4], weight=0.4)

        Use (node, attrdict) tuples to update attributes for specific
        nodes.

        >>> G.add_nodes_from([(1,dict(size=11)), (2,{'color':'blue'})])
        >>> G.node[1]['size']
        11
        >>> H = Graph()
        >>> H.add_nodes_from(G.nodes(data=True))
        >>> H.node[1]['size']
        11

        """
        for n in nodes:
            try:
                newnode=n not in self.succ
            except TypeError:
                nn,ndict = n
                if nn not in self.succ:
                    self.succ[nn] = {}
                    self.pred[nn] = {}
                    newdict = attr.copy()
                    newdict.update(ndict)
                    self.node[nn] = newdict
                else:
                    olddict = self.node[nn]
                    olddict.update(attr)
                    olddict.update(ndict)
                continue
            if newnode:
                self.succ[n] = {}
                self.pred[n] = {}
                self.node[n] = attr.copy()
            else:
                self.node[n].update(attr)

    def remove_node(self, n):
        """Remove node n.

        Removes the node n and all adjacent edges.
        Attempting to remove a non-existent node will raise an exception.

        Parameters
        ----------
        n : node
           A node in the graph

        Raises
        -------
        NetworkXError
           If n is not in the graph.

        See Also
        --------
        remove_nodes_from

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2])
        >>> G.edges()
        [(0, 1), (1, 2)]
        >>> G.remove_node(1)
        >>> G.edges()
        []

        """
        try:
            nbrs=self.succ[n]
            del self.node[n]
        except KeyError: # NetworkXError if n not in self
            raise NetworkXError("The node %s is not in the digraph."%(n,))
        for u in nbrs:
            del self.pred[u][n] # remove all edges n-u in digraph
        del self.succ[n]          # remove node from succ
        for u in self.pred[n]:
            del self.succ[u][n] # remove all edges n-u in digraph
        del self.pred[n]          # remove node from pred


    def remove_nodes_from(self, nbunch):
        """Remove multiple nodes.

        Parameters
        ----------
        nodes : iterable container
            A container of nodes (list, dict, set, etc.).  If a node
            in the container is not in the graph it is silently
            ignored.

        See Also
        --------
        remove_node

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2])
        >>> e = G.nodes()
        >>> e
        [0, 1, 2]
        >>> G.remove_nodes_from(e)
        >>> G.nodes()
        []

        """
        for n in nbunch:
            try:
                succs=self.succ[n]
                del self.node[n]
                for u in succs:
                    del self.pred[u][n] # remove all edges n-u in digraph
                del self.succ[n]          # now remove node
                for u in self.pred[n]:
                    del self.succ[u][n] # remove all edges n-u in digraph
                del self.pred[n]          # now remove node
            except KeyError:
                pass # silent failure on remove


    def add_edge(self, u, v, attr_dict=None, **attr):
        """Add an edge between u and v.

        The nodes u and v will be automatically added if they are
        not already in the graph.

        Edge attributes can be specified with keywords or by providing
        a dictionary with key/value pairs.  See examples below.

        Parameters
        ----------
        u,v : nodes
            Nodes can be, for example, strings or numbers.
            Nodes must be hashable (and not None) Python objects.
        attr_dict : dictionary, optional (default= no attributes)
            Dictionary of edge attributes.  Key/value pairs will
            update existing data associated with the edge.
        attr : keyword arguments, optional
            Edge data (or labels or objects) can be assigned using
            keyword arguments.

        See Also
        --------
        add_edges_from : add a collection of edges

        Notes
        -----
        Adding an edge that already exists updates the edge data.

        Many NetworkX algorithms designed for weighted graphs use as
        the edge weight a numerical value assigned to a keyword
        which by default is 'weight'.

        Examples
        --------
        The following all add the edge e=(1,2) to graph G:

        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> e = (1,2)
        >>> G.add_edge(1, 2)           # explicit two-node form
        >>> G.add_edge(*e)             # single edge as tuple of two nodes
        >>> G.add_edges_from( [(1,2)] ) # add edges from iterable container

        Associate data to edges using keywords:

        >>> G.add_edge(1, 2, weight=3)
        >>> G.add_edge(1, 3, weight=7, capacity=15, length=342.7)
        """
        # set up attribute dict
        if attr_dict is None:
            attr_dict=attr
        else:
            try:
                attr_dict.update(attr)
            except AttributeError:
                raise NetworkXError(\
                    "The attr_dict argument must be a dictionary.")
        # add nodes
        if u not in self.succ:
            self.succ[u]={}
            self.pred[u]={}
            self.node[u] = {}
        if v not in self.succ:
            self.succ[v]={}
            self.pred[v]={}
            self.node[v] = {}
        # add the edge
        datadict=self.adj[u].get(v,{})
        datadict.update(attr_dict)
        self.succ[u][v]=datadict
        self.pred[v][u]=datadict

    def add_edges_from(self, ebunch, attr_dict=None, **attr):
        """Add all the edges in ebunch.

        Parameters
        ----------
        ebunch : container of edges
            Each edge given in the container will be added to the
            graph. The edges must be given as as 2-tuples (u,v) or
            3-tuples (u,v,d) where d is a dictionary containing edge
            data.
        attr_dict : dictionary, optional (default= no attributes)
            Dictionary of edge attributes.  Key/value pairs will
            update existing data associated with each edge.
        attr : keyword arguments, optional
            Edge data (or labels or objects) can be assigned using
            keyword arguments.


        See Also
        --------
        add_edge : add a single edge
        add_weighted_edges_from : convenient way to add weighted edges

        Notes
        -----
        Adding the same edge twice has no effect but any edge data
        will be updated when each duplicate edge is added.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edges_from([(0,1),(1,2)]) # using a list of edge tuples
        >>> e = zip(range(0,3),range(1,4))
        >>> G.add_edges_from(e) # Add the path graph 0-1-2-3

        Associate data to edges

        >>> G.add_edges_from([(1,2),(2,3)], weight=3)
        >>> G.add_edges_from([(3,4),(1,4)], label='WN2898')
        """
        # set up attribute dict
        if attr_dict is None:
            attr_dict=attr
        else:
            try:
                attr_dict.update(attr)
            except AttributeError:
                raise NetworkXError(\
                    "The attr_dict argument must be a dict.")
        # process ebunch
        for e in ebunch:
            ne = len(e)
            if ne==3:
                u,v,dd = e
                assert hasattr(dd,"update")
            elif ne==2:
                u,v = e
                dd = {}
            else:
                raise NetworkXError(\
                    "Edge tuple %s must be a 2-tuple or 3-tuple."%(e,))
            if u not in self.succ:
                self.succ[u] = {}
                self.pred[u] = {}
                self.node[u] = {}
            if v not in self.succ:
                self.succ[v] = {}
                self.pred[v] = {}
                self.node[v] = {}
            datadict=self.adj[u].get(v,{})
            datadict.update(attr_dict)
            datadict.update(dd)
            self.succ[u][v] = datadict
            self.pred[v][u] = datadict


    def remove_edge(self, u, v):
        """Remove the edge between u and v.

        Parameters
        ----------
        u,v: nodes
            Remove the edge between nodes u and v.

        Raises
        ------
        NetworkXError
            If there is not an edge between u and v.

        See Also
        --------
        remove_edges_from : remove a collection of edges

        Examples
        --------
        >>> G = Graph()   # or DiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.remove_edge(0,1)
        >>> e = (1,2)
        >>> G.remove_edge(*e) # unpacks e from an edge tuple
        >>> e = (2,3,{'weight':7}) # an edge with attribute data
        >>> G.remove_edge(*e[:2]) # select first part of edge tuple
        """
        try:
            del self.succ[u][v]
            del self.pred[v][u]
        except KeyError:
            raise NetworkXError("The edge %s-%s not in graph."%(u,v))


    def remove_edges_from(self, ebunch):
        """Remove all edges specified in ebunch.

        Parameters
        ----------
        ebunch: list or container of edge tuples
            Each edge given in the list or container will be removed
            from the graph. The edges can be:

                - 2-tuples (u,v) edge between u and v.
                - 3-tuples (u,v,k) where k is ignored.

        See Also
        --------
        remove_edge : remove a single edge

        Notes
        -----
        Will fail silently if an edge in ebunch is not in the graph.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> ebunch=[(1,2),(2,3)]
        >>> G.remove_edges_from(ebunch)
        """
        for e in ebunch:
            (u,v)=e[:2]  # ignore edge data
            if u in self.succ and v in self.succ[u]:
                del self.succ[u][v]
                del self.pred[v][u]


    def has_successor(self, u, v):
        """Return True if node u has successor v.

        This is true if graph has the edge u->v.
        """
        return (u in self.succ and v in self.succ[u])

    def has_predecessor(self, u, v):
        """Return True if node u has predecessor v.

        This is true if graph has the edge u<-v.
        """
        return (u in self.pred and v in self.pred[u])

    def successors_iter(self,n):
        """Return an iterator over successor nodes of n.

        neighbors_iter() and successors_iter() are the same.
        """
        try:
            return iter(self.succ[n])
        except KeyError:
            raise NetworkXError("The node %s is not in the digraph."%(n,))

    def predecessors_iter(self,n):
        """Return an iterator over predecessor nodes of n."""
        try:
            return iter(self.pred[n])
        except KeyError:
            raise NetworkXError("The node %s is not in the digraph."%(n,))

    def successors(self, n):
        """Return a list of successor nodes of n.

        neighbors() and successors() are the same function.
        """
        return list(self.successors_iter(n))

    def predecessors(self, n):
        """Return a list of predecessor nodes of n."""
        return list(self.predecessors_iter(n))


    # digraph definitions
    neighbors = successors
    neighbors_iter = successors_iter

    def edges_iter(self, nbunch=None, data=False):
        """Return an iterator over the edges.

        Edges are returned as tuples with optional data
        in the order (node, neighbor, data).

        Parameters
        ----------
        nbunch : iterable container, optional (default= all nodes)
            A container of nodes.  The container will be iterated
            through once.
        data : bool, optional (default=False)
            If True, return edge attribute dict in 3-tuple (u,v,data).

        Returns
        -------
        edge_iter : iterator
            An iterator of (u,v) or (u,v,d) tuples of edges.

        See Also
        --------
        edges : return a list of edges

        Notes
        -----
        Nodes in nbunch that are not in the graph will be (quietly) ignored.

        Examples
        --------
        >>> G = DiGraph()   # or MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> [e for e in G.edges_iter()]
        [(0, 1), (1, 2), (2, 3)]
        >>> list(G.edges_iter(data=True)) # default data is {} (empty dict)
        [(0, 1, {}), (1, 2, {}), (2, 3, {})]
        >>> list(G.edges_iter([0,2]))
        [(0, 1), (2, 3)]
        >>> list(G.edges_iter(0))
        [(0, 1)]

        """
        if nbunch is None:
            nodes_nbrs=self.adj.items()
        else:
            nodes_nbrs=((n,self.adj[n]) for n in self.nbunch_iter(nbunch))
        if data:
            for n,nbrs in nodes_nbrs:
                for nbr,data in nbrs.items():
                    yield (n,nbr,data)
        else:
            for n,nbrs in nodes_nbrs:
                for nbr in nbrs:
                    yield (n,nbr)

    # alias out_edges to edges
    out_edges_iter=edges_iter
    out_edges=Graph.edges

    def in_edges_iter(self, nbunch=None, data=False):
        """Return an iterator over the incoming edges.

        Parameters
        ----------
        nbunch : iterable container, optional (default= all nodes)
            A container of nodes.  The container will be iterated
            through once.
        data : bool, optional (default=False)
            If True, return edge attribute dict in 3-tuple (u,v,data).

        Returns
        -------
        in_edge_iter : iterator
            An iterator of (u,v) or (u,v,d) tuples of incoming edges.

        See Also
        --------
        edges_iter : return an iterator of edges
        """
        if nbunch is None:
            nodes_nbrs=self.pred.items()
        else:
            nodes_nbrs=((n,self.pred[n]) for n in self.nbunch_iter(nbunch))
        if data:
            for n,nbrs in nodes_nbrs:
                for nbr,data in nbrs.items():
                    yield (nbr,n,data)
        else:
            for n,nbrs in nodes_nbrs:
                for nbr in nbrs:
                    yield (nbr,n)

    def in_edges(self, nbunch=None, data=False):
        """Return a list of the incoming edges.

        See Also
        --------
        edges : return a list of edges
        """
        return list(self.in_edges_iter(nbunch, data))

    def degree_iter(self, nbunch=None, weight=None):
        """Return an iterator for (node, degree).

        The node degree is the number of edges adjacent to the node.

        Parameters
        ----------
        nbunch : iterable container, optional (default=all nodes)
            A container of nodes.  The container will be iterated
            through once.

        weight : string or None, optional (default=None)
           The edge attribute that holds the numerical value used 
           as a weight.  If None, then each edge has weight 1.
           The degree is the sum of the edge weights adjacent to the node.

        Returns
        -------
        nd_iter : an iterator
            The iterator returns two-tuples of (node, degree).

        See Also
        --------
        degree, in_degree, out_degree, in_degree_iter, out_degree_iter

        Examples
        --------
        >>> G = DiGraph()   # or MultiDiGraph
        >>> G.add_path([0,1,2,3])
        >>> list(G.degree_iter(0)) # node 0 with degree 1
        [(0, 1)]
        >>> list(G.degree_iter([0,1]))
        [(0, 1), (1, 2)]

        """
        if nbunch is None:
            nodes_nbrs=zip(iter(self.succ.items()),iter(self.pred.items()))
        else:
            nodes_nbrs=zip(
                ((n,self.succ[n]) for n in self.nbunch_iter(nbunch)),
                ((n,self.pred[n]) for n in self.nbunch_iter(nbunch)))

        if weight is None:
            for (n,succ),(_,pred) in nodes_nbrs:
                yield (n,len(succ)+len(pred))
        else:
        # edge weighted graph - degree is sum of edge weights
            for (n,succ),(_,pred) in nodes_nbrs:
                yield (n,
                      sum((succ[nbr].get(weight,1) for nbr in succ))+
                      sum((pred[nbr].get(weight,1) for nbr in pred)))


    def in_degree_iter(self, nbunch=None, weight=None):
        """Return an iterator for (node, in-degree).

        The node in-degree is the number of edges pointing in to the node.

        Parameters
        ----------
        nbunch : iterable container, optional (default=all nodes)
            A container of nodes.  The container will be iterated
            through once.

        weight : string or None, optional (default=None)
           The edge attribute that holds the numerical value used 
           as a weight.  If None, then each edge has weight 1.
           The degree is the sum of the edge weights adjacent to the node.

        Returns
        -------
        nd_iter : an iterator
            The iterator returns two-tuples of (node, in-degree).

        See Also
        --------
        degree, in_degree, out_degree, out_degree_iter

        Examples
        --------
        >>> G = DiGraph()
        >>> G.add_path([0,1,2,3])
        >>> list(G.in_degree_iter(0)) # node 0 with degree 0
        [(0, 0)]
        >>> list(G.in_degree_iter([0,1]))
        [(0, 0), (1, 1)]

        """
        if nbunch is None:
            nodes_nbrs=self.pred.items()
        else:
            nodes_nbrs=((n,self.pred[n]) for n in self.nbunch_iter(nbunch))

        if weight is None:
            for n,nbrs in nodes_nbrs:
                yield (n,len(nbrs))
        else:
        # edge weighted graph - degree is sum of edge weights
            for n,nbrs in nodes_nbrs:
                yield (n, sum(data.get(weight,1) for data in nbrs.values()))


    def out_degree_iter(self, nbunch=None, weight=None):
        """Return an iterator for (node, out-degree).

        The node out-degree is the number of edges pointing out of the node.

        Parameters
        ----------
        nbunch : iterable container, optional (default=all nodes)
            A container of nodes.  The container will be iterated
            through once.

        weight : string or None, optional (default=None)
           The edge attribute that holds the numerical value used 
           as a weight.  If None, then each edge has weight 1.
           The degree is the sum of the edge weights adjacent to the node.

        Returns
        -------
        nd_iter : an iterator
            The iterator returns two-tuples of (node, out-degree).

        See Also
        --------
        degree, in_degree, out_degree, in_degree_iter

        Examples
        --------
        >>> G = DiGraph()
        >>> G.add_path([0,1,2,3])
        >>> list(G.out_degree_iter(0)) # node 0 with degree 1
        [(0, 1)]
        >>> list(G.out_degree_iter([0,1]))
        [(0, 1), (1, 1)]

        """
        if nbunch is None:
            nodes_nbrs=self.succ.items()
        else:
            nodes_nbrs=((n,self.succ[n]) for n in self.nbunch_iter(nbunch))

        if weight is None:
            for n,nbrs in nodes_nbrs:
                yield (n,len(nbrs))
        else:
        # edge weighted graph - degree is sum of edge weights
            for n,nbrs in nodes_nbrs:
                yield (n, sum(data.get(weight,1) for data in nbrs.values()))


    def in_degree(self, nbunch=None, weight=None):
        """Return the in-degree of a node or nodes.

        The node in-degree is the number of edges pointing in to the node.

        Parameters
        ----------
        nbunch : iterable container, optional (default=all nodes)
            A container of nodes.  The container will be iterated
            through once.

        weight : string or None, optional (default=None)
           The edge attribute that holds the numerical value used 
           as a weight.  If None, then each edge has weight 1.
           The degree is the sum of the edge weights adjacent to the node.

        Returns
        -------
        nd : dictionary, or number
            A dictionary with nodes as keys and in-degree as values or
            a number if a single node is specified.

        See Also
        --------
        degree, out_degree, in_degree_iter

        Examples
        --------
        >>> G = DiGraph()   # or MultiDiGraph
        >>> G.add_path([0,1,2,3])
        >>> G.in_degree(0)
        0
        >>> G.in_degree([0,1])
        {0: 0, 1: 1}
        >>> list(G.in_degree([0,1]).values())
        [0, 1]
        """
        if nbunch in self:      # return a single node
            return next(self.in_degree_iter(nbunch,weight))[1]
        else:           # return a dict
            return dict(self.in_degree_iter(nbunch,weight))

    def out_degree(self, nbunch=None, weight=None):
        """Return the out-degree of a node or nodes.

        The node out-degree is the number of edges pointing out of the node.

        Parameters
        ----------
        nbunch : iterable container, optional (default=all nodes)
            A container of nodes.  The container will be iterated
            through once.

        weight : string or None, optional (default=None)
           The edge attribute that holds the numerical value used 
           as a weight.  If None, then each edge has weight 1.
           The degree is the sum of the edge weights adjacent to the node.

        Returns
        -------
        nd : dictionary, or number
            A dictionary with nodes as keys and out-degree as values or
            a number if a single node is specified.

        Examples
        --------
        >>> G = DiGraph()   # or MultiDiGraph
        >>> G.add_path([0,1,2,3])
        >>> G.out_degree(0)
        1
        >>> G.out_degree([0,1])
        {0: 1, 1: 1}
        >>> list(G.out_degree([0,1]).values())
        [1, 1]


        """
        if nbunch in self:      # return a single node
            return next(self.out_degree_iter(nbunch,weight))[1]
        else:           # return a dict
            return dict(self.out_degree_iter(nbunch,weight))

    def clear(self):
        """Remove all nodes and edges from the graph.

        This also removes the name, and all graph, node, and edge attributes.

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> G.clear()
        >>> G.nodes()
        []
        >>> G.edges()
        []

        """
        self.succ.clear()
        self.pred.clear()
        self.node.clear()
        self.graph.clear()


    def is_multigraph(self):
        """Return True if graph is a multigraph, False otherwise."""
        return False


    def is_directed(self):
        """Return True if graph is directed, False otherwise."""
        return True

    def to_directed(self):
        """Return a directed copy of the graph.

        Returns
        -------
        G : DiGraph
            A deepcopy of the graph.

        Notes
        -----
        This returns a "deepcopy" of the edge, node, and
        graph attributes which attempts to completely copy
        all of the data and references.

        This is in contrast to the similar D=DiGraph(G) which returns a
        shallow copy of the data.

        See the Python copy module for more information on shallow
        and deep copies, http://docs.python.org/library/copy.html.

        Examples
        --------
        >>> G = Graph()   # or MultiGraph, etc
        >>> G.add_path([0,1])
        >>> H = G.to_directed()
        >>> H.edges()
        [(0, 1), (1, 0)]

        If already directed, return a (deep) copy

        >>> G = DiGraph()   # or MultiDiGraph, etc
        >>> G.add_path([0,1])
        >>> H = G.to_directed()
        >>> H.edges()
        [(0, 1)]
        """
        return deepcopy(self)

    def to_undirected(self, reciprocal=False):
        """Return an undirected representation of the digraph.

        Parameters
        ----------
        reciprocal : bool (optional)
          If True only keep edges that appear in both directions 
          in the original digraph. 

        Returns
        -------
        G : Graph
            An undirected graph with the same name and nodes and
            with edge (u,v,data) if either (u,v,data) or (v,u,data)
            is in the digraph.  If both edges exist in digraph and
            their edge data is different, only one edge is created
            with an arbitrary choice of which edge data to use.
            You must check and correct for this manually if desired.

        Notes
        -----
        If edges in both directions (u,v) and (v,u) exist in the
        graph, attributes for the new undirected edge will be a combination of
        the attributes of the directed edges.  The edge data is updated
        in the (arbitrary) order that the edges are encountered.  For
        more customized control of the edge attributes use add_edge().

        This returns a "deepcopy" of the edge, node, and
        graph attributes which attempts to completely copy
        all of the data and references.

        This is in contrast to the similar G=DiGraph(D) which returns a
        shallow copy of the data.

        See the Python copy module for more information on shallow
        and deep copies, http://docs.python.org/library/copy.html.
        """
        H=Graph()
        H.name=self.name
        H.add_nodes_from(self)
        if reciprocal is True:
            H.add_edges_from( (u,v,deepcopy(d))
                              for u,nbrs in self.adjacency_iter()
                              for v,d in nbrs.items() 
                              if v in self.pred[u])
        else:
            H.add_edges_from( (u,v,deepcopy(d))
                              for u,nbrs in self.adjacency_iter()
                              for v,d in nbrs.items() )
        H.graph=deepcopy(self.graph)
        H.node=deepcopy(self.node)
        return H


    def reverse(self, copy=True):
        """Return the reverse of the graph.

        The reverse is a graph with the same nodes and edges
        but with the directions of the edges reversed.

        Parameters
        ----------
        copy : bool optional (default=True)
            If True, return a new DiGraph holding the reversed edges.
            If False, reverse the reverse graph is created using
            the original graph (this changes the original graph).
        """
        if copy:
            H = self.__class__(name="Reverse of (%s)"%self.name)
            H.add_nodes_from(self)
            H.add_edges_from( (v,u,deepcopy(d)) for u,v,d 
                              in self.edges(data=True) )
            H.graph=deepcopy(self.graph)
            H.node=deepcopy(self.node)
        else:
            self.pred,self.succ=self.succ,self.pred
            self.adj=self.succ
            H=self
        return H


    def subgraph(self, nbunch):
        """Return the subgraph induced on nodes in nbunch.

        The induced subgraph of the graph contains the nodes in nbunch
        and the edges between those nodes.

        Parameters
        ----------
        nbunch : list, iterable
            A container of nodes which will be iterated through once.

        Returns
        -------
        G : Graph
            A subgraph of the graph with the same edge attributes.

        Notes
        -----
        The graph, edge or node attributes just point to the original graph.
        So changes to the node or edge structure will not be reflected in
        the original graph while changes to the attributes will.

        To create a subgraph with its own copy of the edge/node attributes use:
        Graph(G.subgraph(nbunch))

        If edge attributes are containers, a deep copy can be obtained using:
        G.subgraph(nbunch).copy()

        For an inplace reduction of a graph to a subgraph you can remove nodes:
        G.remove_nodes_from([ n in G if n not in set(nbunch)])

        Examples
        --------
        >>> G = Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_path([0,1,2,3])
        >>> H = G.subgraph([0,1,2])
        >>> H.edges()
        [(0, 1), (1, 2)]
        """
        bunch = self.nbunch_iter(nbunch)
        # create new graph and copy subgraph into it
        H = self.__class__()
        # namespace shortcuts for speed
        H_succ=H.succ
        H_pred=H.pred
        self_succ=self.succ
        # add nodes
        for n in bunch:
            H_succ[n]={}
            H_pred[n]={}
        # add edges
        for u in H_succ:
            Hnbrs=H_succ[u]
            for v,datadict in self_succ[u].items():
                if v in H_succ:
                    # add both representations of edge: u-v and v-u
                    Hnbrs[v]=datadict
                    H_pred[v][u]=datadict
        # copy node and attribute dictionaries
        for n in H:
            H.node[n]=self.node[n]
        H.graph=self.graph
        return H

'''
============================================================
convert module, from networkx
============================================================
'''
"""
This module provides functions to convert 
NetworkX graphs to and from other formats.

The preferred way of converting data to a NetworkX graph 
is through the graph constuctor.  The constructor calls
the to_networkx_graph() function which attempts to guess the
input type and convert it automatically.

Examples
--------

Create a 10 node random graph from a numpy matrix

>>> import numpy
>>> a=numpy.reshape(numpy.random.random_integers(0,1,size=100),(10,10))
>>> D=DiGraph(a) 

or equivalently

>>> D=to_networkx_graph(a,create_using=DiGraph()) 

Create a graph with a single edge from a dictionary of dictionaries

>>> d={0: {1: 1}} # dict-of-dicts single edge (0,1)
>>> G=Graph(d)


See Also
--------
nx_pygraphviz, nx_pydot

"""
#    Copyright (C) 2006-2011 by 
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

__all__ = ['to_networkx_graph',
           'from_dict_of_dicts', 'to_dict_of_dicts',
           'from_dict_of_lists', 'to_dict_of_lists',
           'from_edgelist', 'to_edgelist',
           'from_numpy_matrix', 'to_numpy_matrix',
           'to_numpy_recarray',
           'from_scipy_sparse_matrix', 'to_scipy_sparse_matrix']

def _prep_create_using(create_using):
    """Return a graph object ready to be populated.

    If create_using is None return the default (just networkx.Graph())
    If create_using.clear() works, assume it returns a graph object.
    Otherwise raise an exception because create_using is not a networkx graph.

    """
    if create_using is None:
        G=Graph()
    else:
        G=create_using
        try:
            G.clear()
        except:
            raise TypeError("Input graph is not a networkx graph type")
    return G

def to_networkx_graph(data,create_using=None,multigraph_input=False):
    """Make a NetworkX graph from a known data structure.

    The preferred way to call this is automatically
    from the class constructor

    >>> d={0: {1: {'weight':1}}} # dict-of-dicts single edge (0,1)
    >>> G=Graph(d)
    
    instead of the equivalent

    >>> G=from_dict_of_dicts(d)

    Parameters
    ----------
    data : a object to be converted
       Current known types are:
         any NetworkX graph
         dict-of-dicts
         dist-of-lists
         list of edges
         numpy matrix
         numpy ndarray
         scipy sparse matrix
         pygraphviz agraph

    create_using : NetworkX graph
       Use specified graph for result.  Otherwise a new graph is created.

    multigraph_input : bool (default False)
      If True and  data is a dict_of_dicts,
      try to create a multigraph assuming dict_of_dict_of_lists.
      If data and create_using are both multigraphs then create
      a multigraph from a multigraph.

    """
    # NX graph
    if hasattr(data,"adj"):
        try:
            result= from_dict_of_dicts(data.adj,\
                    create_using=create_using,\
                    multigraph_input=data.is_multigraph())
            if hasattr(data,'graph') and isinstance(data.graph,dict):
                result.graph=data.graph.copy()
            if hasattr(data,'node') and isinstance(data.node,dict):
                result.node=dict( (n,dd.copy()) for n,dd in data.node.items() )
            return result
        except:
            raise NetworkXError("Input is not a correct NetworkX graph.")

    # pygraphviz  agraph
#    if hasattr(data,"is_strict"):
#        try:
#            return from_agraph(data,create_using=create_using)
#        except:
#            raise NetworkXError("Input is not a correct pygraphviz graph.")

    # dict of dicts/lists
    if isinstance(data,dict):
        try:
            return from_dict_of_dicts(data,create_using=create_using,\
                    multigraph_input=multigraph_input)
        except:
            try:
                return from_dict_of_lists(data,create_using=create_using)
            except:
                raise TypeError("Input is not known type.")

    # list or generator of edges
    if (isinstance(data,list)
        or hasattr(data,'next')
        or hasattr(data, '__next__')): 
        try:
            return from_edgelist(data,create_using=create_using)
        except:
            raise NetworkXError("Input is not a valid edge list")

    raise NetworkXError(\
          "Input is not a known data type for conversion.")

    return 


def convert_to_undirected(G):
    """Return a new undirected representation of the graph G.

    """
    return G.to_undirected()


def convert_to_directed(G):
    """Return a new directed representation of the graph G.

    """
    return G.to_directed()


def to_dict_of_lists(G,nodelist=None):
    """Return adjacency representation of graph as a dictionary of lists.

    Parameters
    ----------
    G : graph
       A NetworkX graph 

    nodelist : list       
       Use only nodes specified in nodelist

    Notes
    -----
    Completely ignores edge data for MultiGraph and MultiDiGraph.

    """
    if nodelist is None:
        nodelist=G

    d = {}
    for n in nodelist:
        d[n]=[nbr for nbr in G.neighbors(n) if nbr in nodelist]
    return d            

def from_dict_of_lists(d,create_using=None):
    """Return a graph from a dictionary of lists.

    Parameters
    ----------
    d : dictionary of lists
      A dictionary of lists adjacency representation.

    create_using : NetworkX graph
       Use specified graph for result.  Otherwise a new graph is created.

    Examples
    --------
    >>> dol= {0:[1]} # single edge (0,1)
    >>> G=from_dict_of_lists(dol)

    or
    >>> G=Graph(dol) # use Graph constructor

    """
    G=_prep_create_using(create_using)
    G.add_nodes_from(d)        
    if G.is_multigraph() and not G.is_directed():
        # a dict_of_lists can't show multiedges.  BUT for undirected graphs,
        # each edge shows up twice in the dict_of_lists.  
        # So we need to treat this case separately.
        seen={}
        for node,nbrlist in d.items():
            for nbr in nbrlist:
                if nbr not in seen:
                    G.add_edge(node,nbr)
            seen[node]=1  # don't allow reverse edge to show up 
    else:
        G.add_edges_from( ((node,nbr) for node,nbrlist in d.items() 
                           for nbr in nbrlist) )
    return G                         


def to_dict_of_dicts(G,nodelist=None,edge_data=None):
    """Return adjacency representation of graph as a dictionary of dictionaries.

    Parameters
    ----------
    G : graph
       A NetworkX graph 

    nodelist : list       
       Use only nodes specified in nodelist

    edge_data : list, optional       
       If provided,  the value of the dictionary will be
       set to edge_data for all edges.  This is useful to make
       an adjacency matrix type representation with 1 as the edge data.
       If edgedata is None, the edgedata in G is used to fill the values.
       If G is a multigraph, the edgedata is a dict for each pair (u,v).
    
    """
    dod={}
    if nodelist is None:
        if edge_data is None:
            for u,nbrdict in G.adjacency_iter():
                dod[u]=nbrdict.copy()
        else: # edge_data is not None
            for u,nbrdict in G.adjacency_iter():
                dod[u]=dod.fromkeys(nbrdict, edge_data)
    else: # nodelist is not None
        if edge_data is None:
            for u in nodelist:
                dod[u]={}
                for v,data in ((v,data) for v,data in G[u].items() if v in nodelist):
                    dod[u][v]=data
        else: # nodelist and edge_data are not None
            for u in nodelist:
                dod[u]={}
                for v in ( v for v in G[u] if v in nodelist):
                    dod[u][v]=edge_data
    return dod

def from_dict_of_dicts(d,create_using=None,multigraph_input=False):
    """Return a graph from a dictionary of dictionaries.

    Parameters
    ----------
    d : dictionary of dictionaries
      A dictionary of dictionaries adjacency representation.

    create_using : NetworkX graph
       Use specified graph for result.  Otherwise a new graph is created.

    multigraph_input : bool (default False)
       When True, the values of the inner dict are assumed 
       to be containers of edge data for multiple edges.
       Otherwise this routine assumes the edge data are singletons.

    Examples
    --------
    >>> dod= {0: {1:{'weight':1}}} # single edge (0,1)
    >>> G=from_dict_of_dicts(dod)

    or
    >>> G=Graph(dod) # use Graph constructor

    """
    G=_prep_create_using(create_using)
    G.add_nodes_from(d)
    # is dict a MultiGraph or MultiDiGraph?
    if multigraph_input:
        # make a copy of the list of edge data (but not the edge data)
        if G.is_directed():  
            if G.is_multigraph():
                G.add_edges_from( (u,v,key,data)
                                  for u,nbrs in d.items() 
                                  for v,datadict in nbrs.items() 
                                  for key,data in datadict.items()
                                )
            else:
                G.add_edges_from( (u,v,data)
                                  for u,nbrs in d.items() 
                                  for v,datadict in nbrs.items() 
                                  for key,data in datadict.items()
                                )
        else: # Undirected
            if G.is_multigraph():
                seen=set()   # don't add both directions of undirected graph
                for u,nbrs in d.items():
                    for v,datadict in nbrs.items():
                        if (u,v) not in seen:
                            G.add_edges_from( (u,v,key,data) 
                                               for key,data in datadict.items()
                                              )
                            seen.add((v,u)) 
            else:
                seen=set()   # don't add both directions of undirected graph
                for u,nbrs in d.items():
                    for v,datadict in nbrs.items():
                        if (u,v) not in seen:
                            G.add_edges_from( (u,v,data)
                                        for key,data in datadict.items() )
                            seen.add((v,u)) 

    else: # not a multigraph to multigraph transfer
        if G.is_multigraph() and not G.is_directed():
            # d can have both representations u-v, v-u in dict.  Only add one.
            # We don't need this check for digraphs since we add both directions,
            # or for Graph() since it is done implicitly (parallel edges not allowed)
            seen=set()
            for u,nbrs in d.items():
                for v,data in nbrs.items():
                    if (u,v) not in seen:
                        G.add_edge(u,v,attr_dict=data)
                    seen.add((v,u))
        else:
            G.add_edges_from( ( (u,v,data) 
                                for u,nbrs in d.items() 
                                for v,data in nbrs.items()) )
    return G                         

def to_edgelist(G,nodelist=None):
    """Return a list of edges in the graph.

    Parameters
    ----------
    G : graph
       A NetworkX graph 

    nodelist : list       
       Use only nodes specified in nodelist

    """
    if nodelist is None:
        return G.edges(data=True)
    else:
        return G.edges(nodelist,data=True)

def from_edgelist(edgelist,create_using=None):
    """Return a graph from a list of edges.

    Parameters
    ----------
    edgelist : list or iterator
      Edge tuples 

    create_using : NetworkX graph
       Use specified graph for result.  Otherwise a new graph is created.

    Examples
    --------
    >>> edgelist= [(0,1)] # single edge (0,1)
    >>> G=from_edgelist(edgelist)

    or
    >>> G=Graph(edgelist) # use Graph constructor

    """
    G=_prep_create_using(create_using)
    G.add_edges_from(edgelist)
    return G                         

def to_numpy_matrix(G, nodelist=None, dtype=None, order=None,
                    multigraph_weight=sum, weight='weight'):
    """Return the graph adjacency matrix as a NumPy matrix.

    Parameters
    ----------
    G : graph
        The NetworkX graph used to construct the NumPy matrix.

    nodelist : list, optional       
       The rows and columns are ordered according to the nodes in `nodelist`.
       If `nodelist` is None, then the ordering is produced by G.nodes().

    dtype : NumPy data type, optional
        A valid single NumPy data type used to initialize the array. 
        This must be a simple type such as int or numpy.float64 and
        not a compound data type (see to_numpy_recarray)
        If None, then the NumPy default is used.

    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory. If None, then the NumPy default 
        is used.

    multigraph_weight : {sum, min, max}, optional        
        An operator that determines how weights in multigraphs are handled.
        The default is to sum the weights of the multiple edges.

    weight : string or None   optional (default='weight')
        The edge attribute that holds the numerical value used for 
        the edge weight.  If None then all edge weights are 1.


    Returns
    -------
    M : NumPy matrix
       Graph adjacency matrix.

    See Also
    --------
    to_numpy_recarray, from_numpy_matrix

    Notes
    -----
    The matrix entries are assigned with weight edge attribute. When
    an edge does not have the weight attribute, the value of the entry is 1.
    For multiple edges, the values of the entries are the sums of the edge
    attributes for each edge.

    When `nodelist` does not contain every node in `G`, the matrix is built 
    from the subgraph of `G` that is induced by the nodes in `nodelist`.
    
    Examples
    --------
    >>> G = MultiDiGraph()
    >>> G.add_edge(0,1,weight=2)
    >>> G.add_edge(1,0)
    >>> G.add_edge(2,2,weight=3)
    >>> G.add_edge(2,2)
    >>> to_numpy_matrix(G, nodelist=[0,1,2])
    matrix([[ 0.,  2.,  0.],
            [ 1.,  0.,  0.],
            [ 0.,  0.,  4.]])

    """
    try:
        import numpy as np
    except ImportError:
        raise ImportError(\
          "to_numpy_matrix() requires numpy: http://scipy.org/ ")

    if nodelist is None:
        nodelist = G.nodes()

    nodeset = set(nodelist)
    if len(nodelist) != len(nodeset):
        msg = "Ambiguous ordering: `nodelist` contained duplicates."
        raise NetworkXError(msg)

    nlen=len(nodelist)
    undirected = not G.is_directed()
    index=dict(zip(nodelist,range(nlen)))

    if G.is_multigraph():
        # Handle MultiGraphs and MultiDiGraphs
        # array of nan' to start with, any leftover nans will be converted to 0
        # nans are used so we can use sum, min, max for multigraphs
        M = np.zeros((nlen,nlen), dtype=dtype, order=order)+np.nan
        # use numpy nan-aware operations
        operator={sum:np.nansum, min:np.nanmin, max:np.nanmax}
        try:
            op=operator[multigraph_weight]
        except:
            raise ValueError('multigraph_weight must be sum, min, or max')

        for u,v,attrs in G.edges_iter(data=True):
            if (u in nodeset) and (v in nodeset):
                i,j = index[u],index[v]
                e_weight = attrs.get(weight, 1)
                M[i,j] = op([e_weight,M[i,j]]) 
                if undirected:
                    M[j,i] = M[i,j]
        # convert any nans to zeros
        M = np.asmatrix(np.nan_to_num(M))
    else:
        # Graph or DiGraph, this is much faster than above 
        M = np.zeros((nlen,nlen), dtype=dtype, order=order)
        for u,nbrdict in G.adjacency_iter():
            for v,d in nbrdict.items():
                try:
                    M[index[u],index[v]]=d.get(weight,1)
                except KeyError:
                    pass
        M = np.asmatrix(M)
    return M


def from_numpy_matrix(A,create_using=None):
    """Return a graph from numpy matrix.

    The numpy matrix is interpreted as an adjacency matrix for the graph.

    Parameters
    ----------
    A : numpy matrix
      An adjacency matrix representation of a graph

    create_using : NetworkX graph
       Use specified graph for result.  The default is Graph()

    Notes
    -----
    If the numpy matrix has a single data type for each matrix entry it 
    will be converted to an appropriate Python data type.  

    If the numpy matrix has a user-specified compound data type the names
    of the data fields will be used as attribute keys in the resulting 
    NetworkX graph.

    See Also
    --------
    to_numpy_matrix, to_numpy_recarray

    Examples
    --------
    Simple integer weights on edges:

    >>> import numpy
    >>> A=numpy.matrix([[1,1],[2,1]])
    >>> G=from_numpy_matrix(A)

    User defined compound data type on edges:

    >>> import numpy
    >>> dt=[('weight',float),('cost',int)]
    >>> A=numpy.matrix([[(1.0,2)]],dtype=dt)                      
    >>> G=from_numpy_matrix(A)
    >>> G.edges(data=True)
    [(0, 0, {'cost': 2, 'weight': 1.0})]
    """
    kind_to_python_type={'f':float,
                         'i':int,
                         'u':int,
                         'b':bool,
                         'c':complex,
                         'S':str,
                         'V':'void'}

    try: # Python 3.x
        _ = chr(1245) # just to trigger the exception
        kind_to_python_type['U']=str
    except ValueError: # Python 2.6+
        kind_to_python_type['U']=unicode

    # This should never fail if you have created a numpy matrix with numpy...  
    try:
        import numpy as np
    except ImportError:
        raise ImportError(\
          "from_numpy_matrix() requires numpy: http://scipy.org/ ")

    G=_prep_create_using(create_using)
    n,m=A.shape
    if n!=m:
        raise NetworkXError("Adjacency matrix is not square.",
                               "nx,ny=%s"%(A.shape,))
    dt=A.dtype
    try:
        python_type=kind_to_python_type[dt.kind]
    except:
        raise TypeError("Unknown numpy data type: %s"%dt)

    # make sure we get isolated nodes
    G.add_nodes_from(range(n)) 
    # get a list of edges
    x,y=np.asarray(A).nonzero()         

    # handle numpy constructed data type
    if python_type is 'void':
        fields=sorted([(offset,dtype,name) for name,(dtype,offset) in
                       A.dtype.fields.items()])
        for (u,v) in zip(x,y):         
            attr={}
            for (offset,dtype,name),val in zip(fields,A[u,v]):
                attr[name]=kind_to_python_type[dtype.kind](val)
            G.add_edge(u,v,attr)
    else: # basic data type
        G.add_edges_from( ((u,v,{'weight':python_type(A[u,v])}) 
                           for (u,v) in zip(x,y)) )
    return G


def to_numpy_recarray(G,nodelist=None,
                      dtype=[('weight',float)],
                      order=None):
    """Return the graph adjacency matrix as a NumPy recarray.

    Parameters
    ----------
    G : graph
        The NetworkX graph used to construct the NumPy matrix.

    nodelist : list, optional       
       The rows and columns are ordered according to the nodes in `nodelist`.
       If `nodelist` is None, then the ordering is produced by G.nodes().

    dtype : NumPy data-type, optional
        A valid NumPy named dtype used to initialize the NumPy recarray. 
        The data type names are assumed to be keys in the graph edge attribute 
        dictionary.

    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory. If None, then the NumPy default 
        is used.

    Returns
    -------
    M : NumPy recarray
       The graph with specified edge data as a Numpy recarray 

    Notes
    -----
    When `nodelist` does not contain every node in `G`, the matrix is built 
    from the subgraph of `G` that is induced by the nodes in `nodelist`.
    
    Examples
    --------
    >>> G = Graph()
    >>> G.add_edge(1,2,weight=7.0,cost=5)
    >>> A=to_numpy_recarray(G,dtype=[('weight',float),('cost',int)])
    >>> print(A.weight)
    [[ 0.  7.]
     [ 7.  0.]]
    >>> print(A.cost)
    [[0 5]
     [5 0]]
    """
    try:
        import numpy as np
    except ImportError:
        raise ImportError(\
          "to_numpy_matrix() requires numpy: http://scipy.org/ ")

    if G.is_multigraph():
        raise NetworkXError("Not implemented for multigraphs.")

    if nodelist is None:
        nodelist = G.nodes()

    nodeset = set(nodelist)
    if len(nodelist) != len(nodeset):
        msg = "Ambiguous ordering: `nodelist` contained duplicates."
        raise NetworkXError(msg)

    nlen=len(nodelist)
    undirected = not G.is_directed()
    index=dict(zip(nodelist,range(nlen)))
    M = np.zeros((nlen,nlen), dtype=dtype, order=order)

    names=M.dtype.names
    for u,v,attrs in G.edges_iter(data=True):
        if (u in nodeset) and (v in nodeset):
            i,j = index[u],index[v]
            values=tuple([attrs[n] for n in names])
            M[i,j] = values
            if undirected:
                M[j,i] = M[i,j]

    return M.view(np.recarray)


def to_scipy_sparse_matrix(G, nodelist=None, dtype=None, 
                           weight='weight', fmt='csr'):
    """Return the graph adjacency matrix as a SciPy sparse matrix.

    Parameters
    ----------
    G : graph
        The NetworkX graph used to construct the NumPy matrix.

    nodelist : list, optional       
       The rows and columns are ordered according to the nodes in `nodelist`.
       If `nodelist` is None, then the ordering is produced by G.nodes().

    dtype : NumPy data-type, optional
        A valid NumPy dtype used to initialize the array. If None, then the
        NumPy default is used.

    weight : string or None   optional (default='weight')
        The edge attribute that holds the numerical value used for 
        the edge weight.  If None then all edge weights are 1.

    fmt : str in {'bsr', 'csr', 'csc', 'coo', 'lil', 'dia', 'dok'} 
        The type of the matrix to be returned (default 'csr').  For
        some algorithms different implementations of sparse matrices
        can perform better.  See [1]_ for details.
    
    Returns
    -------
    M : SciPy sparse matrix
       Graph adjacency matrix.

    Notes
    -----
    The matrix entries are populated using the edge attribute held in 
    parameter weight. When an edge does not have that attribute, the 
    value of the entry is 1.

    For multiple edges the matrix values are the sums of the edge weights.

    When `nodelist` does not contain every node in `G`, the matrix is built 
    from the subgraph of `G` that is induced by the nodes in `nodelist`.
    
    Uses lil_matrix format. To convert to other formats specify the 
    format= keyword.

    Examples
    --------
    >>> G = MultiDiGraph()
    >>> G.add_edge(0,1,weight=2)
    >>> G.add_edge(1,0)
    >>> G.add_edge(2,2,weight=3)
    >>> G.add_edge(2,2)
    >>> S = to_scipy_sparse_matrix(G, nodelist=[0,1,2])
    >>> S.todense()
    matrix([[ 0.,  2.,  0.],
            [ 1.,  0.,  0.],
            [ 0.,  0.,  4.]])
    
    References
    ----------
    .. [1] Scipy Dev. References, 
       "Sparse Matrices"  
       http://docs.scipy.org/doc/scipy/reference/sparse.html
    """
    try:
        from scipy import sparse
    except ImportError:
        raise ImportError(\
          "to_scipy_sparse_matrix() requires scipy: http://scipy.org/ ")

    if nodelist is None:
        nodelist = G.nodes()

    nodeset = set(nodelist)
    if len(nodelist) != len(nodeset):
        msg = "Ambiguous ordering: `nodelist` contained duplicates."
        raise NetworkXError(msg)

    nlen=len(nodelist)
    undirected = not G.is_directed()
    index=dict(zip(nodelist,range(nlen)))
    M = sparse.lil_matrix((nlen,nlen), dtype=dtype)

    for u,v,attrs in G.edges_iter(data=True):
        if (u in nodeset) and (v in nodeset):
            i,j = index[u],index[v]
            M[i,j] += attrs.get(weight, 1)
            if undirected:
                M[j,i] = M[i,j]
    try:
        return M.asformat(fmt)
    except AttributeError:
        raise NetworkXError("Unknown sparse matrix format: %s"%fmt)
    

def from_scipy_sparse_matrix(A,create_using=None):
    """Return a graph from scipy sparse matrix adjacency list. 

    Parameters
    ----------
    A : scipy sparse matrix
      An adjacency matrix representation of a graph

    create_using : NetworkX graph
       Use specified graph for result.  The default is Graph()

    Examples
    --------
    >>> import scipy.sparse
    >>> A=scipy.sparse.eye(2,2,1)
    >>> G=from_scipy_sparse_matrix(A)

    """
    G=_prep_create_using(create_using)

    # convert all formats to lil - not the most efficient way       
    AA=A.tolil()
    n,m=AA.shape

    if n!=m:
        raise NetworkXError(\
              "Adjacency matrix is not square. nx,ny=%s"%(A.shape,))
    G.add_nodes_from(range(n)) # make sure we get isolated nodes

    for i,row in enumerate(AA.rows):
        for pos,j in enumerate(row):
            G.add_edge(i,j,**{'weight':AA.data[i][pos]})
    return G

'''
============================================================
Algorithms for directed acyclic graphs (DAGs), from networkx
============================================================
'''

# Exception handling

# the root of all Exceptions
class NetworkXException(Exception):
    """Base class for exceptions in NetworkX."""

class NetworkXError(NetworkXException):
    """Exception for a serious error in NetworkX"""

class NetworkXAlgorithmError(NetworkXException):
    """Exception for unexpected termination of algorithms."""

class NetworkXUnfeasible(NetworkXAlgorithmError):
    """Exception raised by algorithms trying to solve a problem
    instance that has no feasible solution."""

"""Algorithms for directed acyclic graphs (DAGs)."""
#    Copyright (C) 2006-2011 by 
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
__author__ = """\n""".join(['Aric Hagberg <aric.hagberg@gmail.com>',
                            'Dan Schult (dschult@colgate.edu)',
                            'Ben Edwards (bedwards@cs.unm.edu)'])
__all__ = ['topological_sort', 
           'topological_sort_recursive',
           'is_directed_acyclic_graph',
           'is_aperiodic']

def is_directed_acyclic_graph(G):
    """Return True if the graph G is a directed acyclic graph (DAG) or 
    False if not.
    
    Parameters
    ----------
    G : NetworkX graph
      A graph

    Returns
    -------
    is_dag : bool
       True if G is a DAG, false otherwise
    """
    try:
        topological_sort(G)
        return True
    except NetworkXUnfeasible:
        return False

def topological_sort(G,nbunch=None):
    """Return a list of nodes in topological sort order.

    A topological sort is a nonunique permutation of the nodes
    such that an edge from u to v implies that u appears before v in the
    topological sort order.

    Parameters
    ----------
    G : NetworkX digraph
       A directed graph

    nbunch : container of nodes (optional)
       Explore graph in specified order given in nbunch

    Raises
    ------
    NetworkXError
       Topological sort is defined for directed graphs only. If the
       graph G is undirected, a NetworkXError is raised.

    NetworkXUnfeasible
       If G is not a directed acyclic graph (DAG) no topological sort
       exists and a NetworkXUnfeasible exception is raised.

    Notes
    -----
    This algorithm is based on a description and proof in
    The Algorithm Design Manual [1]_ .

    See also
    --------
    is_directed_acyclic_graph

    References
    ----------
    .. [1] Skiena, S. S. The Algorithm Design Manual  (Springer-Verlag, 1998). 
        http://www.amazon.com/exec/obidos/ASIN/0387948600/ref=ase_thealgorithmrepo/
    """
    if not G.is_directed():
        raise NetworkXError(
                "Topological sort not defined on undirected graphs.")

    # nonrecursive version
    seen={}
    order_explored=[] # provide order and 
    explored={}       # fast search without more general priorityDictionary
                     
    if nbunch is None:
        nbunch = G.nodes_iter() 
    for v in nbunch:     # process all vertices in G
        if v in explored: 
            continue
        fringe=[v]   # nodes yet to look at
        while fringe:
            w=fringe[-1]  # depth first search
            if w in explored: # already looked down this branch
                fringe.pop()
                continue
            seen[w]=1     # mark as seen
            # Check successors for cycles and for new nodes
            new_nodes=[]
            for n in G[w]:
                if n not in explored:
                    if n in seen: #CYCLE !!
                        raise NetworkXUnfeasible("Graph contains a cycle.")
                    new_nodes.append(n)
            if new_nodes:   # Add new_nodes to fringe
                fringe.extend(new_nodes)
            else:           # No new nodes so w is fully explored
                explored[w]=1
                order_explored.insert(0,w) # reverse order explored
                fringe.pop()    # done considering this node
    return order_explored

def topological_sort_recursive(G,nbunch=None):
    """Return a list of nodes in topological sort order.

    A topological sort is a nonunique permutation of the nodes such
    that an edge from u to v implies that u appears before v in the
    topological sort order.

    Parameters
    ----------
    G : NetworkX digraph

    nbunch : container of nodes (optional)
       Explore graph in specified order given in nbunch

    Raises
    ------
    NetworkXError
       Topological sort is defined for directed graphs only. If the
       graph G is undirected, a NetworkXError is raised.

    NetworkXUnfeasible
        If G is not a directed acyclic graph (DAG) no topological sort
        exists and a NetworkXUnfeasible exception is raised.

    Notes
    -----
    This is a recursive version of topological sort.

    See also
    --------
    topological_sort
    is_directed_acyclic_graph

    """
    if not G.is_directed():
        raise NetworkXError(
                "Topological sort not defined on undirected graphs.")

    # function for recursive dfs
    def _dfs(G,seen,explored,v):
        seen.add(v)
        for w in G[v]:
            if w not in seen: 
                if not _dfs(G,seen,explored,w):
                    return False
            elif w in seen and w not in explored:
                # cycle Found--- no topological sort
                raise NetworkXUnfeasible("Graph contains a cycle.")
        explored.insert(0,v) # inverse order of when explored 
        return True

    seen=set()
    explored=[]

    if nbunch is None:
        nbunch = G.nodes_iter() 
    for v in nbunch:  # process all nodes
        if v not in explored:
            if not _dfs(G,seen,explored,v): 
                raise NetworkXUnfeasible("Graph contains a cycle.")
    return explored

def is_aperiodic(G):
    """Return True if G is aperiodic.

    A directed graph is aperiodic if there is no integer k > 1 that 
    divides the length of every cycle in the graph.

    Parameters
    ----------
    G : NetworkX DiGraph
      Graph

    Returns
    -------
    aperiodic : boolean
      True if the graph is aperiodic False otherwise

    Raises
    ------
    NetworkXError
      If G is not directed

    Notes
    -----
    This uses the method outlined in [1]_, which runs in O(m) time
    given m edges in G. Note that a graph is not aperiodic if it is
    acyclic as every integer trivial divides length 0 cycles.

    References
    ----------
    .. [1] Jarvis, J. P.; Shier, D. R. (1996),
       Graph-theoretic analysis of finite Markov chains,
       in Shier, D. R.; Wallenius, K. T., Applied Mathematical Modeling:
       A Multidisciplinary Approach, CRC Press.
    """
    if not G.is_directed():
        raise NetworkXError("is_aperiodic not defined for undirected graphs")

    s = next(G.nodes_iter())
    levels = {s:0}
    this_level = [s]
    g = 0
    l = 1
    while this_level:
        next_level = []
        for u in this_level:
            for v in G[u]:
                if v in levels: # Non-Tree Edge
                    g = gcd(g, levels[u]-levels[v] + 1)
                else: # Tree Edge
                    next_level.append(v)
                    levels[v] = l
        this_level = next_level
        l += 1
    if len(levels)==len(G): #All nodes in tree
        return g==1
    else:
        return g==1 and is_aperiodic(G.subgraph(set(G)-set(levels)))

"""
=================================
Depth-first search, from networkx 
=================================

Basic algorithms for depth-first searching.

Based on http://www.ics.uci.edu/~eppstein/PADS/DFS.py
by D. Eppstein, July 2004.
"""
__author__ = """\n""".join(['Aric Hagberg <hagberg@lanl.gov>'])

__all__ = ['dfs_edges', 'dfs_tree',
           'dfs_predecessors', 'dfs_successors',
           'dfs_preorder_nodes','dfs_postorder_nodes',
           'dfs_labeled_edges']

def dfs_edges(G,source=None):
    """Produce edges in a depth-first-search starting at source."""
    # Based on http://www.ics.uci.edu/~eppstein/PADS/DFS.py
    # by D. Eppstein, July 2004.
    if source is None:
        # produce edges for all components
        nodes=G
    else:
        # produce edges for components with source
        nodes=[source]
    visited=set()
    for start in nodes:
        if start in visited:
            continue
        visited.add(start)
        stack = [(start,iter(G[start]))]
        while stack:
            parent,children = stack[-1]
            try:
                child = next(children)
                if child not in visited:
                    yield parent,child
                    visited.add(child)
                    stack.append((child,iter(G[child])))
            except StopIteration:
                stack.pop()


def dfs_tree(G, source=None):
    """Return directed tree of depth-first-search from source."""
    return DiGraph(dfs_edges(G,source=source))


def dfs_predecessors(G, source=None):
    """Return dictionary of predecessors in depth-first-search from source."""
    return dict((t,s) for s,t in dfs_edges(G,source=source))


def dfs_successors(G, source=None):
    """Return dictionary of successors in depth-first-search from source."""
    d=defaultdict(list)
    for s,t in dfs_edges(G,source=source):
        d[s].append(t)
    return dict(d)


def dfs_postorder_nodes(G,source=None):
    """Produce nodes in a depth-first-search post-ordering starting 
    from source.
    """
    post=(v for u,v,d in dfs_labeled_edges(G,source=source)
          if d['dir']=='reverse')
    # chain source to end of pre-ordering
#    return chain(post,[source])
    return post


def dfs_preorder_nodes(G,source=None):
    """Produce nodes in a depth-first-search pre-ordering starting at source."""
    pre=(v for u,v,d in dfs_labeled_edges(G,source=source) 
         if d['dir']=='forward')
    # chain source to beginning of pre-ordering
#    return chain([source],pre)
    return pre


def dfs_labeled_edges(G,source=None):
    """Produce edges in a depth-first-search starting at source and
    labeled by direction type (forward, reverse, nontree).
    """
    # Based on http://www.ics.uci.edu/~eppstein/PADS/DFS.py
    # by D. Eppstein, July 2004.
    if source is None:
        # produce edges for all components
        nodes=G
    else:
        # produce edges for components with source
        nodes=[source]
    visited=set()
    for start in nodes:
        if start in visited:
            continue
        yield start,start,{'dir':'forward'}
        visited.add(start)
        stack = [(start,iter(G[start]))]
        while stack:
            parent,children = stack[-1]
            try:
                child = next(children)
                if child in visited:
                    yield parent,child,{'dir':'nontree'}
                else:
                    yield parent,child,{'dir':'forward'}
                    visited.add(child)
                    stack.append((child,iter(G[child])))
            except StopIteration:
                stack.pop()
                if stack:
                    yield stack[-1][0],parent,{'dir':'reverse'}
        yield start,start,{'dir':'reverse'}

'''
============================================================
Pedtplot pedigree data structure.

Both individuals and marriage nodes are treated as Person
objects, which are the nodes of the pedigree's graph member.
============================================================
'''
class Pedigree(object):
    def __init__(self, graph, data):
        '''Construct a pedigree from a Person DAG (edge goes from parent to child).'''
        
        # Validate graph
        if (not is_directed_acyclic_graph(graph)):
            raise ValueError('Pedigree graph must be a directed acyclic graph (DAG)')
        for node in graph.nodes_iter():
            degree = graph.in_degree(node)
            if (degree > 2):
                raise ValueError('In-degree must be at most 2, but degree[%d] = %d' % (node, degree, ))
        self.data = data
        
        # Convert pedigree to the form we need: add dummy founder and marriage nodes. Keep
        # two graphs: person_graph has persons only; graph has persons and marriages
        self.person_graph = self.__add_dummy_founders(graph)
        self.graph = self.__add_marriage_nodes(self.person_graph)
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return 'Pedigree[nodes=%d, edges=%d]' % (self.graph.number_of_nodes(), self.graph.number_of_edges()) 

    def pprint(self):
        print 'Pedigree[nodes=%d, edges=%d]' % (self.graph.number_of_nodes(), self.graph.number_of_edges())
        print 'Node F  M  Sex  Type  Label'
        for (name, person) in self.data.iteritems():
            print '%6d %6d %6d %4d %4d %s' % (name, person.father, person.mother, 
                                           person.sex, person.node_type, person.label) 

    def father(self, person):
        '''Father node.'''
        return self.data[person.father] if person.father != MISSING else None

    def mother(self, person):
        '''Mother node.'''
        return self.data[person.mother] if person.mother != MISSING else None
   
    def nodes_with_type(self, node_type):
        '''Return all nodes of type node_type.''' 
        return (k for (k, v) in self.data.iteritems() if v.node_type == node_type)

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def undirected_graph(self):
        '''Return the undirected pedigree graph.''' 
        return self.graph.to_undirected()

    @property
    def person_nodes(self):
        '''Return all person nodes.'''
        return self.person_graph.nodes_iter()

    @property
    def marriage_nodes(self):
        '''Return all marriage nodes.''' 
        return self.nodes_with_type(NODE_TYPE.MARRIAGE)

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __add_dummy_founders(self, g):
        '''Add dummy founders to children with a single parent. The pedfiddler algorithm only supports 0
        or 2 parents.'''
        g_new = g.copy()
        
        # Increment dummy node names starting after the largest ID in the current pedigree
        dummy_name = max(name for name in g.nodes_iter())
        for name in g.nodes_iter():
            person = self.data[name]
            f, m = person.father, person.mother
            if ((f == MISSING) ^ (m == MISSING)):
                # If node has a single parent, create dummy node for the other parent
                dummy_name += 1
                dummy_parent = Person(dummy_name, node_type=NODE_TYPE.DUMMY,
                                      sex=SEX.MALE if (f == MISSING) else SEX.FEMALE)
                g_new.add_node(dummy_name)
                g_new.add_edge(dummy_name, name)
                self.data[dummy_name] = dummy_parent
                if f == MISSING: person.father = dummy_name                    
                else: person.mother = dummy_name
        return g_new

    def __add_marriage_nodes(self, g):
        '''Return an extended pedigree graph that includes marriage nodes and edges.
        Marriage node IDs are negative to avoid collision with real nodes.'''
        g_new =  g.copy()
        name = -1 #  
        for (parents, children) in Pedigree.__families(g):
            father = [x for x in parents if self.data[x].sex == SEX.MALE][0] 
            mother = [x for x in parents if self.data[x].sex == SEX.FEMALE][0]
            marriage = Person(name, node_type=NODE_TYPE.MARRIAGE, father=father, mother=mother)
            self.data[name] = marriage
            g_new.remove_edges_from((p, c) for p in parents for c in children) 
            g_new.add_edges_from((parent, name) for parent in parents)
            g_new.add_edges_from((name, child) for child in children)
            name -= 1
        return g_new

    @staticmethod
    def __families(g):
        '''A generator of all marriages: tuples (parents, children) in a DAG g.'''
        visited = set()
        # Loop over nodes, considering each as a child. Child -> parents -> all children that
        # share the same parents -> output this as a family -> mark children as visited. This
        # treats polygamous pedigrees as well.
        for child in g.nodes():
            if not child in visited:
                parents = set(g.predecessors(child))
                # Child belongs to a family only if it is not a founder
                if parents:
                    children = set(node for node in all_successors(g, parents)
                                   if set(g.predecessors(node)) == parents)
                    visited |= children
                    yield (parents, children)

#========================================================================================
class Person(object):
    '''An immutable object representing a haploid organism (e.g., a person). Serves as the
    node type of a Pedigree graph.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, name, node_type, sex=SEX.UNKNOWN, status=STATUS.UNKNOWN,
                 father=MISSING, mother=MISSING, label=None):
        '''Construct a person with a unique identifier 'name' and descriptive fields.'''
        # Descriptive attributes
        self.name       = name
        self.node_type  = node_type
        self.sex        = sex
        self.status     = status
        self.father     = father
        self.mother     = mother
        self.label      = label if label else name  

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return 'Person[%d]' % (self.name,) 

    def __key(self):
        '''Hash key.'''
        return self.name

    def __eq__(self, other):
        '''Equality of objects.'''
        return self.__key() == other.__key()

    def __ne__(self, other):
        '''Inequality of objects.'''
        return self.__key() != other.__key()

    def __cmp__(self, other):
        return cmp(self.name, other.name)
    
    def copy(self):
        '''Return a deep copy of this object.'''
        return Person(name=self.name, node_type=self.node_type, sex=self.sex, status=self.status,
                       father=self.father, mother=self.mother, label=self.label)

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def parents(self):
        '''Return the parent IDs of this node.''' 
        return (self.father, self.mother)

#========================================================================================
class PedigreeInfo(object):
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, ped, data=None):
        '''Construct an empty pedigree information object. Requires the pedigree graph.'''
        self.ped = ped
        self.data = data if data else dict((name, PersonInfo(name)) for name in ped.graph.nodes_iter())
   
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return self.pprint()

    def pprint(self, max_gen=None):
        g = self.ped.graph
        s  = 'PedigreeInfo[nodes=%d, edges=%d]\n' % (g.number_of_nodes(), g.number_of_edges())
        s += 'Node     F      M      Sex   Type    Label    Gen    Coord    Color\n'
        for (name, person) in self.ped.data.iteritems():
            info = self.data[name]
            if max_gen is None or info.gen <= max_gen:
                s += '%6d %6d %6d %4d %4d %2d (%.2f,%.2f) (%3d,%3d,%3d)\n' % \
                ((name, person.father, person.mother, person.sex, person.node_type) +
                 (info.gen,) + info.coord + info.rgb)
        return s

    def to_txt(self, max_gen=None):
        s = ''
        for (name, person) in self.ped.data.iteritems():
            info = self.data[name]
            if (max_gen is None or info.gen <= max_gen) and person.node_type != NODE_TYPE.MARRIAGE:
                s += '%6d %6d %6d %4d %4d %-7s %2d %.2f,%.2f %3d %3d %3d\n' % \
                ((name, person.father, person.mother, person.sex, person.node_type, person.label) +
                 (info.gen,) + info.coord + info.rgb)
        return s

    def copy(self):
        '''Return a deep copy of this object (graph is shallow-copied, data is deep-copied).'''
        return PedigreeInfo(self.ped, self.data.copy())
  
    def marriage_data_at_level(self, gen):
        '''Return the list of marriage nodes just below generation gen.'''
        (g, data) = (self.graph, self.data)
        return [data[m] for m in self.marriage_nodes
                if max(data[x].gen for x in g.predecessors(m)) == gen]

    def person_data_at_level(self, gen):
        '''Return the list of person nodes at generation gen.'''
        data = self.data
        return [data[x] for x in self.person_nodes if data[x].gen == gen]

    def coord(self, index=None):
        '''Vector of coordinates of all nodes in graph node iteration order.'''
        return [self.data[name].coord[index] for name in self.ped.graph.nodes_iter()] if index is not None else [self.data[name].coord for name in self.graph.nodes_iter()]

    def set_coord(self, coord):
        '''Set coordinate data from a dictionary of node-to-positions.'''
        for (name, value) in coord.iteritems():
            self.data[name].coord = value

    def set_rgb(self, colors):
        '''Set node colors from a dictionary of affection-status-to-color.'''
        (node, info) = (self.ped.data, self.data)
        for name in self.person_nodes:
            info[name].rgb  = colors[node[name].status]
    
    def balance_marriages(self):
        for node in self.ped.marriage_nodes:
            self.__balance_marriage(self.ped.data[node], self.data[node])
    
    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def graph(self):
        '''Return the full (marriage) graph.'''
        return self.ped.graph

    @property
    def person_graph(self):
        '''Return the person graph.'''
        return self.ped.person_graph

    @property
    def person_nodes(self):
        '''Return all person nodes.'''
        return self.ped.person_nodes

    @property
    def marriage_nodes(self):
        '''Return all marriage nodes.''' 
        return self.ped.marriage_nodes

    @property
    def gen_dict(self):
        '''Dictionary of node generation numbers.'''
        return dict((name, self.data[name].gen) for name in self.ped.person_graph.nodes_iter())

    @property
    def coord_dict(self):
        '''Dictionary of node coordinates.'''
        return dict((name, self.data[name].coord) for name in self.ped.graph.nodes_iter())

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __balance_marriage(self, node, info):
        '''Set marriage node's x-coordinate to be on the line connecting the center of mass
        of the parents and the children.'''
        (f, m) = (self.data[node.father], self.data[node.mother])
        x1 = 0.5*(f.coord[0] + m.coord[0])
        y1 = 0.5*(f.coord[1] + m.coord[1])
        
        children = self.ped.graph.successors(node.name)
        x2 = average(self.data[child].coord[0] for child in children)
        y2 = average(self.data[child].coord[1] for child in children)

        yold = info.coord[1]
        info.coord = (x1 - (y1 - yold)/(y1 - y2)*(x1 - x2), yold)

#========================================================================================
class PersonInfo(object):
    '''A mutable object representing data attributes about a Person.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, name, gen=0, degree=0, coord=(-1.0, -1.0), rgb=(255,255,255)):
        '''Construct a person with a unique identifier name [and fields].'''
        # Descriptive attributes
        self.name       = name
        self.gen        = gen
        self.degree     = degree
        self.coord      = coord
        self.rgb        = rgb
        #self.flag = const.DOWN
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return 'PersonInfo[name=%d, coord=(%d,%d), rgb=(%d,%d,%d)]' % \
            ((self.name,) + self.coord + self.rgb)
            
    def copy(self):
        '''Return a deep copy of this object.'''
        return PersonInfo(name=self.name, gen=self.gen, degree=self.degree,
                          coord=self.coord, rgb=self.rgb)

########################################################################
# Pedigree I/O
########################################################################
def read_pedigree(infile, separator=' '):
    '''Read pedigree from text file.'''
    # Build nodes and a dictionary of id-to-node
    data = dict((person.name, person) for person in __read_person_lines(infile, separator))
    g = DiGraph()
    g.add_nodes_from(data.iterkeys())
    
    # Build edges; add a dummy entry to node dictionary for a missing node ID in edges
    g.add_edges_from(((parent, child) for (child, parent) in
                      reduce(tuple.__add__, (((name, person.father), (name, person.mother)) 
                                             for (name, person) in data.iteritems()))
                     if parent != MISSING))
    return Pedigree(g, data)

def __read_person_lines(infile, separator):
    '''Read pedigree file row-by-row and emit Node objects.'''
        # Read delimited data from standard input
    with open(infile, 'rb') as infile:
        reader = csv.reader(infile, delimiter=separator, skipinitialspace=True)
        # Extract requested columns from relevant lines
        for line in reader:
            if line: # Make sure there's at least one entry
                # Expected row format: <id> <father_id> <mother_id> <sex> <affection_status> <label>
                yield Person(name=int(line[0]),
                             node_type=NODE_TYPE.PERSON,
                             father=int(line[1]),
                             mother=int(line[2]),
                             sex=int(line[3]), 
                             status=int(line[4]),
                             label=line[5])

'''
#===================================================================
inline istream & operator >>(istream &s, rawpedrecord &p)
{
  char buff[1025];
  char kid[128];
  char fa[128];
  char mo[128];
  int i;
  int x1, x2;
  s.getline(buff,1024);    
  istringstream in(buff);
  in >> kid >> fa >> mo;
  if(_sex_flag) in >> x1;
  else x1=0;
  if(_status_flag) in >> x2;
  else x2=0;
  if(_label_flag){
    i = in.tellg();
    if(i == -1){
      i=0; buff[0] = 0;
    }
    while(buff[i] == ' ' || buff[i] == '\t') i++;
  }else{
    for(i = 0; buff[i] != ' '  && buff[i] != '\t' && buff[i] != 0; i++);
    buff[i] = 0; i = 0;
  }  
  p.update(kid,fa,mo,x1,x2,&buff[i]);
  return s;
}


void rawpedrecord::update(char *Id, char *Dad, char *Mom, int _sex, int _status,char *lab){
    strcpy(id,Id);
    if(Dad != NULL){
       dad = new char[strlen(Dad)+1];
    strcpy(dad,Dad);
    }
    if(Mom != NULL){
       mom = new char[strlen(Mom)+1];
    strcpy(mom,Mom);
    }
    sex = _sex; status = _status;
    if(lab != NULL){
      strcpy(label,lab);
    }
}

void check1(const int n, const rawpedrecord *p){
  for(int i = 0; i<n; i++){
    if(strcmp(p[i].id,p[i].dad) == 0 || strcmp(p[i].id, p[i].mom) == 0){
      if(_console_)
    cerr << "ERROR: '" << p[i].id << "' cannot be its own parent." <<  endl;
      else
    fl_alert("ERROR: '%s' cannot be its own parent.",p[i].id);
      sys.exit(1);
    }
  }
  for(int n1 = n-1, i=0; i<n1; i++)
    for(int j=i+1; j<n; j++)
      if(strcmp(p[i].id, p[j].id) == 0){
    if(_console_)
      cerr << "ERROR: '" << p[i].id << "' is not unique." <<  endl;
    else
      fl_alert("ERROR: '%s' is not unique.", p[i].id);
    sys.exit(1);
      }
} 

void FixPrecedence(const int n, rawpedrecord *p)
{
  check1(n,p);
  rawpedrecord tmp;
  int i,j, n_1, jo;
  n_1 = n-1;
  for(i=0; i <n_1 ; i++){
    jo = 2*(n-i);
    for(j=i+1; j < n; j++){
      if(strcmp(p[i].dad, p[j].id) == 0){
    if(!jo--){
      if(_console_)
        cerr << "ERROR: '" << p[j].id << "' cannot be its own ancestor."<<endl;
      else
        fl_alert("ERROR: '%s' cannot be its own ancestor.",p[j].id);
      sys.exit(1);
    }
    tmp = p[j];
    p[j] = p[i];
    p[i] = tmp;
    j = i;
    continue;
      }
      if(strcmp(p[i].mom, p[j].id) == 0){
    if(!jo--){
      if(_console_)
        cerr << "ERROR: '" <<p[j].id << "' cannot be its own ancestor." <<endl;
      else
        fl_alert("ERROR: '%s' cannot be its own ancestor.",p[j].id);
      sys.exit(1);
    }       
    tmp = p[j];
    p[j] = p[i];
    p[i] = tmp;
    j = i;
      }
    }
  }
}

bool isnotfounder(const char *parent){
  if(parent){
    const char *p = parent;
    while(*p == ' ') p++;
    if(p[0] == '\0' || 
       (p[0] == '0' && p[1] == '\0')) return 0;
  }
  return 1;
}

void pedigree::append(const char *Id, const char *Fa, const char *Mo,
              const int &_sex, const int &_status, const char *lab){
  unsigned int i,j;
  if((i = isnotfounder(Fa))){
    for(; i<= N; i++) 
      if(strcmp(Fa,id[i]) == 0) break;
    if(i > N)
      add(Fa,0,0,1,0,"");
  }
  if((j = isnotfounder(Mo))){
    for(; j<= N; j++) 
      if(strcmp(Mo,id[j]) == 0) break;
    if(j > N)
      add(Mo,0,0,2,0,"");
  }
  add(Id, i, j, _sex, _status,lab);
}

pedigree::pedigree(const char *filename){
    FixPrecedence(n,rec);
    setup( n + 100 );
    for(i=0; i< n; i++)
      append( rec[i].id, rec[i].dad, rec[i].mom, rec[i].sex, rec[i].status,
              rec[i].label );
  }
}

void readlist_ (const char *filename,int sexes, int affection,int labels,
        int rn, int gn, int bn, int rd, int gd, int bd,
        int ru, int gu, int bu)
{
  int rows, cols;
  int ierror = 1;
  _sex_flag =sexes;
  _status_flag=affection;
  _label_flag=labels;
  unsigned char r,g,b;
  rows = numlines(filename, cols);
  pedigree _ped(filename);

  initpedc(2*_ped.N+1,_ped.N+1,2*_ped.N+1,_ped.N+1,2*_ped.N+1);
  
  for(unsigned int j=1; j <=_ped.N; j++){
    switch(_ped.status[j]){
    case 2:
      r=rd; g=gd; b=bd;
      break;
    case 0:
      r=ru; g=gu; b=bu;
      break;
    default:
      r=rn; g=gn; b=bn;
    }
    ierror &= setindc(j,_ped.father[j],_ped.mother[j],
              _ped.sex[j],r,g,b,_ped.id[j],_ped.label[j]);
  }
  if(!ierror){ 
    if(_console_)
      cerr << "Pedigree indexing error." << endl;
    else
      fl_alert("Pedigree indexing error.");
    sys.exit(1);
  }
  pedset=1;
  for (int i=1; i<=indpt; i++) indstack[i].info.flag = UP;
  for (int i=1; i<=marpt; i++) marstack[i].info.flag = UP;
  mngrcoords();
  coordset =1;
}

void initpedc(int mxind,int mxmar,int mxindlist,int mxmarlist, int mxlabel)

   Allocate stacks used to store the pedigree.
   
   mxind     = length of the array of individuals
   mxmar     = length of the array of marriages
   mxindlist = length of the array of indlists
   mxmarlist = length of the array of marlists
   mxlabel   = length of the ped array = largest number used as an identifier

{
int i;

maxind = mxind;
maxmar = mxmar;
maxindlist = mxindlist;
maxmarlist = mxmarlist;
maxlabel = mxlabel;

indstack = (struct ind *)calloc((unsigned)maxind+1,sizeof(*indstack));

indliststack = (struct indlist *)calloc((unsigned)maxindlist+1,sizeof(*indliststack));
marstack = (struct mar *)calloc((unsigned)maxmar+1,sizeof(*marstack));
marliststack = (struct marlist *)calloc((unsigned)maxmarlist+1,sizeof(*marliststack));

ped = (struct ind **)calloc((unsigned)maxlabel+1,sizeof(*ped));
ped[0] = indstack;
for(i=1; i < maxlabel; i++) ped[i]=NULL;

nomar = marstack;
marstack[0].pa = indstack;
marstack[0].ma = indstack;
persons[0] = Person(None)
}

struct inddat {
  int name, sex, gen, degree;
  int flag, index;
  float xcoord, ycoord;
  unsigned char r,g,b;
  char *id;
  char *label;
};

struct ind { 
  struct ind *pa, *ma; 
  struct mar *par; 
  struct marlist *own; 
  struct inddat info;
};

struct mardat{
  int flag;
  float xcoord, ycoord;
};

struct mar{
  struct ind *pa, *ma;
  struct indlist *kids;
  struct mardat info;
};

#ifdef __cplusplus
extern "C" {
#endif

/*==================================================*/
/* External calls. See pedstacks.c for definitions. */
/*==================================================*/

  extern int maxind, maxmar, maxindlist, maxmarlist, maxlabel,
    indpt, marpt, indlistpt, marlistpt,
    higen, hilabel, ncomp, 
    coordset, genreset, pedset;
  
  extern struct ind **ped, *indstack, *noind;
  extern struct mar *marstack, *nomar;
  noindlist = [];
  nomarlist = [];
'''
                
def all_successors(g, nodeset):
    '''Return the set of sucessor nodes of the node set nodeset in the DAG g.'''
    result = set([])
    for node in nodeset:
        result |= set(g.successors(node))
    return result

def average(a):
    '''Return the average value of the collection a.'''
    (s, count) = (0.0, 0)
    for x in a:
        s += x
        count += 1
    return s/(1.0*count) if count > 0 else None

'''
============================================================
Node coordinate computation.

Corresponds to pedfiddler files:
annealcoord.c
============================================================
'''

########################################################################
# Constants
########################################################################

'''Default node colors based on affection status.'''
DEFAULT_COLORS = {STATUS.UNKNOWN : (255,255,255),
                  STATUS.NORMAL  : (255,  0,  0), 
                  STATUS.AFFECTED: (255,222,173)}

########################################################################
# Public Methods
########################################################################
class CoordParams(object):
    def __init__(self):
        self.algorithm          = 'hardcoded' # 'annealing'
        self.balance_marriages  = True  # Balance marriages in post-processing
        
        # simulated annealing parameters
        self.num_iterations     = 1000
        self.initial_temp       = 15.0
        self.cooling_factor     = 0.9
        
        # Hard-coded coordinates
        self.coord_file         = None
    
def compute_coords(ped_info, options):
    '''Main method to switch between coordinate computation algorithms.'''
    if options.algorithm == 'hardcoded':
        __CoordComputerHardcoded().run(ped_info, options.coord_file)
    elif options.algorithm == 'default':
        __CoordComputerDefault(ped_info).run(options)
    else:
        raise ValueError('Unsupported coordinate computation algorithm ''%s''' % (options.algorithm,))

    if options.balance_marriages:
        # Optional: seems to place marriages better instead of on top of each other
        ped_info.balance_marriages()

#========================================================================================
class __CoordComputerHardcoded():
    '''Read node coordinates from a delimited file and set them in the ped_info object.'''
    ########################################################################
    # Public Methods
    ########################################################################            
    def run(self, ped_info, in_file):
        '''Compute coordinates and set them in ped_info.'''
        ped_info.set_coord(dict((k,v) for (k,v) in self.__read_coord_lines(in_file)))

    ########################################################################
    # Private Methods
    ########################################################################                
    def __read_coord_lines(self, infile):
        '''Read node coordinates from a delimited file and set them in the ped_info object.'''
        with open(infile, 'rb') as f:
            reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
            for line in reader:
                if line:
                    yield (int(line[0]), (float(line[1]), float(line[2])))

#========================================================================================
class __CoordComputerDefault():
    ########################################################################
    # Constructors
    ########################################################################            
    def __init__(self, ped_info):
        '''Compute coordinates and set them in ped_info.'''
        self.ped_info = ped_info
        self.node_info = ped_info.data

    ########################################################################
    # Public Methods
    ########################################################################            
    def run(self, options):
        '''Compute coordinates and set them in ped_info.'''
        self.__set_generations()
        #print self.ped_info.gen_dict
        self.__set_coord()

    ########################################################################
    # Private Methods
    ########################################################################                
    def __set_generations(self):
        '''Computes generation numbers so that later children's generation number is one more than
        both parents' generation numbers. Founders at the top of the pedigree are set to generation 0.'''
        (g, node_info) = (self.ped_info.person_graph, self.ped_info.data)
        
        # Initialize generation to 1 at all nodes, and in particular top nodes
        for name in g.nodes_iter():
            node_info[name].gen = 0

        # Since a child is processed after both its parents, set its generation # according to theirs    
        for name in reversed(list(dfs_postorder_nodes(g))):
            parents = g.predecessors(name)
            if parents:
                node_info[name].gen = max(node_info[x].gen for x in parents)+1
            #print name, 'parents', parents, 'gen', node_info[name].gen, 'parents gen', [node_info[x].gen for x in parents]

        # Advance parent's generation so that it is one less than his/her earliest child generation
        for name in dfs_postorder_nodes(g):
            children = g.successors(name)
            if children:
                node_info[name].gen = min(node_info[x].gen for x in children)-1
            #print name, 'children', children, 'gen', node_info[name].gen, 'children gen', [node_info[x].gen for x in children]  
            if node_info[name].gen < 0:
                raise ValueError('Negative gen at node %s' % (name,))

    def __set_coord(self):
        '''Assigns coordinates for a marriage node graph of the pedigree.'''
        num_generations = max(self.node_info[x].gen for x in self.ped_info.person_graph.nodes_iter())
        y_unit = 0.5/(num_generations+1)
        
        '''New order by coordinates: left to right, top to bottom '''
        gen = 0
        marriages = self.__place_top_marriages(gen, 1 - 1.0/(num_generations+1))
        self.__place_top_persons(marriages, gen, y_unit)
 
        for a in xrange(1, num_generations+1):
            persons = self.__place_non_founders(a, y_unit)
            founders = self.__place_founders(a, 1-(a+1.0)/(num_generations+1.0) + y_unit)
            #print 'founders', founders
            
            # Sort generation by x-coordinate (left to right) and replace it at equally spaced points.
            # Main loop is over spouses; place founder next to each spouse if exists.
            sorted_persons = ()
            for x in [i[0] for i in sorted(((info.name, info.coord[0]) for info in persons), key=lambda x: x[1])]:
                if founders.has_key(x):
                    sorted_persons += (x, founders[x])
                else:
                    sorted_persons += (x,)
            #print 'Sorted persons', sorted_persons

            scale_factor = 1.0/(len(sorted_persons)+1)
            for (count, name) in enumerate(sorted_persons):
                info = self.node_info[name]
                info.coord = ((count+1)*scale_factor, info.coord[1]) 
                #print name, info.name, info.coord

            if a < num_generations:
                self.__place_marriages(a, 1-(a+1.0)/(num_generations+1.0))

    def __place_top_marriages(self, gen, y):
        '''Enter the first level marriages into marstack and set their coordinates.'''
        #print '__place_top_marriages', gen
        marriages = self.ped_info.marriage_data_at_level(gen)
        
        for (count, info) in enumerate(marriages):
            info.coord = (count+1.0, y)
            
        scale_factor = 1.0/(len(marriages)+1.0)
        for info in marriages:
            info.coord = (info.coord[0]*scale_factor, info.coord[1]) 
        return marriages
    
    def __place_top_persons(self, marriages, gen, y_unit):
        '''Set the coordinates of generation-0 persons.'''
        #print '__place_top_persons', gen
        (g, node_info) = (self.ped_info.graph, self.node_info)
        #persons = self.ped_info.person_data_at_level(gen)
        persons = []
        for a in ([node_info[x] for x in g.predecessors(info.name)] for info in marriages):
            persons += a
        
        for (count, info) in enumerate(persons):
            '''Compute xcoord as the average of its marriages '''
            name = info.name
            marriages = g.successors(name)
            if not marriages:
                raise ValueError('Gen 0 individual %d has no marriages!' % (name, ))
            else:
                info.coord = (average(node_info[x].coord[0] for x in marriages), 
                              max(node_info[x].coord[1] for x in marriages) + y_unit)
                                
        scale_factor = 1.0/(len(persons) + 1.0)
        for (count, info) in enumerate(persons): 
            info.coord = (scale_factor*(count+1.0), info.coord[1])
            #print info.name, info.coord
    
    def __place_non_founders(self, gen, y_unit):
        '''Place non-founder individuals of generation gen.'''
        #print '__place_non_founders', gen
        (g, node_info) = (self.ped_info.graph, self.node_info)
        persons = [x for x in self.ped_info.person_data_at_level(gen) if len(g.predecessors(x.name)) > 0]
        for info in persons:
            '''Place node directly below its preceding marriage node.'''
            name = info.name
            marriage_coord = node_info[g.predecessors(name)[0]].coord 
            info.coord = (marriage_coord[0], marriage_coord[1] - y_unit)
            #print info.name, info.coord
        return persons
    
    def __place_founders(self, gen, y):
        '''Place founder individuals of generation gen level's y-coordinate y. Return a dictionary
        of spouse-to-founder dciontary list. 
        
        NOTE: At this time, we do not support polygamous pedigrees, and
        each founder is assumed to have a unique marriage and spouse.'''
        #print '__place_founders', gen
        g = self.ped_info.graph
        persons = [x for x in self.ped_info.person_data_at_level(gen) if len(g.predecessors(x.name)) == 0]
        for info in persons:
            '''Place node at the generation gen level. Another option is to use its spouse's y-coord.'''
            info.coord = (info.coord[0], y)
            #print info.name, info.coord
        # This is where we assume there's a unique successor (marriage) and spouse
        founders = {}
        for x in [info.name for info in persons]:
            p = g.predecessors(g.successors(x)[0]) 
            s = [spouse for spouse in p if spouse != x][0]
            founders[s] = x
        return founders
        
    def __place_marriages(self, gen, y):
        '''Place marriages of generation gen. Each marriage's x-coordinate is the average
        of its parents' x-coordinates.'''
        #print '__place_marriages', gen
        (g, node_info) = (self.ped_info.graph, self.node_info)
        for info in self.ped_info.marriage_data_at_level(gen):
            #print g.predecessors(info.name)
            #print [node_info[x].coord[0] for x in g.predecessors(info.name)]
            info.coord = (average(node_info[x].coord[0] for x in g.predecessors(info.name)), y)
            #print info.name, info.coord
  
#========================================================================================
#
## Annealing - not supported yet
#
###include 'headio.h'
###include 'Headped.H'
###include 'headrandom.h'
###include 'util.h'
##
##int *ningen;
##int **marpic;
##
##void annealcoords(int nit, float boil, float cool)
##{
##  ''' 
##   Uses the annealing algorithm to minimise the total horizontal components of 
##     line length for a marriage node graph of the pedigree.
##     The perturbation scheme is to swap the coordinates of two marriages and 
##     adjust the individuals involved accordingly.
##  '''
##
##  {
##    struct mar *m;
##    int i, j, a, nin, oldpt, most;
##    float spacer, y_unit;
##    ningen = (int *)calloc(num_generations+1,sizeof(int));
##    marpic = (int **)calloc(num_generations+1,sizeof(int *));
##    y_unit = .5/num_generations;
##
##    for (i=1; i<=marpt; i++) {
##      a = MAXIMUM(marstack[i].pa.info.gen,marstack[i].ma.info.gen);
##      ningen[a]++;
##    }
##    for (i=1, most=0; i<num_generations; i++){
##      marpic[i] = (int *)calloc(2*(ningen[i]+1),sizeof(int));
##
##      if (most < ningen[i]) most = ningen[i];
##      ningen[i] = 0;
##    }
##    spacer = 1.0/((most+1.0)*8.0);
##    for (i=1; i<=marpt; i++){
##      a = MAXIMUM(marstack[i].pa.info.gen,marstack[i].ma.info.gen);
##      marpic[a][++ningen[a]] = i;
##    }
##
##
##    oldpt=marpt;
##
##    for (i=1; i<num_generations; i++) {
##
##    nin = ningen[i];
##    for (j=nin+1; j<=2*nin; j++) {
##      ningen[i]++;
##      m = nextmar();
##      marpic[i][ningen[i]] = marpt;
##      m.info.ycoord = marstack[marpic[i][1]].info.ycoord;
##      a = marpic[i][ningen[i]-nin];
##      m.info.xcoord = marstack[a].info.xcoord;
##      m.info.xcoord += 1.0/(2*(nin+1.0));
##    }
##      
##    }
##    marpt = oldpt;
##    tryswitching(nit,boil,cool,spacer,y_unit);   
##  }
##}
## 
##'''======================'''
##''' Annealing procedures '''
##'''======================'''
##
##define INDCOORDS(X) \
##{\
##if ((X).own != nomarlist)\
##  {\
##  float tot = 0.0;\
##  int num = 0;\
##  struct marlist *ml;\
##  for (ml = (X).own; ml != nomarlist; ml = ml.next)\
##    {\
##    tot += ml.cur.info.xcoord;\
##    num ++;\
##    }\
##  (X).info.xcoord = tot/num;\
##  if ((X).par != nomar)\
##    {\
##    (X).info.xcoord = ((X).info.xcoord + (X).par.info.xcoord) / 2.0 ;\
##    }\
##  }\
##else\
##  (X).info.xcoord = (X).par.info.xcoord;\
##}
##
##define DISTABOUT(X) \
##{\
##distab = 0.0;\
##if ((X).par != nomar)\
##  {\
##  xt = (X).info.xcoord - (X).par.info.xcoord;\
##  distab += xt*xt ;\
##  }\
##for (ml = (X).own; ml != nomarlist; ml = ml.next)\
##  {\
##  xt = (X).info.xcoord - ml.cur.info.xcoord;\
##  distab += xt*xt ;\
##  }\
##}
##
##define COPYANDDEC(A) \
##{\
##newpt++;\
##indstack[newpt].pa = A;\
##indstack[newpt].info.xcoord = A.info.xcoord;\
##indstack[newpt].info.ycoord = A.info.ycoord;\
##DISTABOUT(A)\
##curtot -= distab;\
##}
##
##define SWITCHCOORDS(A,B) \
##{\
##xt = A.info.xcoord;\
##yt = A.info.ycoord;\
##A.info.xcoord = B.info.xcoord;\
##A.info.ycoord = B.info.ycoord;\
##B.info.xcoord = xt;\
##B.info.ycoord = yt;\
##if (A.pa != noind)\
##  {\
##  if (A.pa.info.degree > 1) INDCOORDS(A.pa)\
##  if (A.ma.info.degree > 1) INDCOORDS(A.ma)\
##  for( il = A.kids; il != noindlist; il = il.next)\
##    if (il.cur.info.degree > 1) INDCOORDS(il.cur)\
##  }\
##if (B.pa != noind)\
##  {\
##  if (B.pa.info.degree > 1) INDCOORDS(B.pa)\
##  if (B.ma.info.degree > 1) INDCOORDS(B.ma)\
##  for( il = B.kids; il != noindlist; il = il.next)\
##    if (il.cur.info.degree > 1) INDCOORDS(il.cur)\
##  }\
##}
##
##'''==================='''
##''' Annealing program '''
##'''==================='''
##void tryswitching(int ntries, float boiling, float cooler, 
##          float spacer, float y_unit)
##{
##
##''' Annealing parameters '''
##int nup = 0, ndown = 0, nacross = 0, lastup = 0, lastdown = 0, lastacross = 0;
## float temp = boiling;
##float freezing=0.0;
##float curtot = 0, oldtot;
##
##
##''' Indexes '''
##int i, j, c, d, e;
##int newpt;
##int ndiscon;
##struct marlist *ml;
##struct indlist *il;
##struct mar *a, *b;
##float xt, yt, distab;
##
##
##''' Calculate initial length '''
##for (i=1; i<=indpt; i++) 
##  if (indstack[i].info.degree > 1) { DISTABOUT(&indstack[i]) curtot += distab; }
##
##''' Output initial data '''
##fprintf(stderr,'\n\n Searching for a minimum total squared line length graph \n');
##fprintf(stderr,'\n Total squared x distances at start       = %12.5f',curtot);
##fprintf(stderr,'\n Starting temperature                     = %12.5f',boiling);
##fprintf(stderr,'\n Cooling factor                           = %12.5f',cooler);
##fflush(stderr);
##fprintf(stderr,'\n ');
##
##''' Main loop '''
## randinit();
##for (i = 1; i <= ntries; i++, temp *= cooler)
##  {
##
##  ''' Select marriages to switch '''
##  do { c = randint(1,num_generations-1); } while (ningen[c] < 2);
##  do
##    {
##    d = randint(1,ningen[c]);
##    do { e = randint(1,ningen[c]); } while ( e == d );
##    d = marpic[c][d];
##    e = marpic[c][e];
##    }
##  while ( d > marpt and e > marpt );
##  a = &marstack[d];
##  b = &marstack[e];
##
##  ''' Switch the marriages '''
##  oldtot = curtot;
##  newpt = indpt;
##  if (a.pa != noind)
##    {
##    if (a.pa.info.degree > 1) COPYANDDEC(a.pa)
##    if (a.ma.info.degree > 1) COPYANDDEC(a.ma)
##      for (il= a.kids; il != noindlist; il = il.next)
##      if (il.cur.info.degree > 1) COPYANDDEC(il.cur)
##      }
##  if (b.pa != noind)
##    {
##    if (b.pa.info.degree > 1) COPYANDDEC(b.pa)
##    if (b.ma.info.degree > 1) COPYANDDEC(b.ma)
##     for (il= b.kids; il != noindlist; il = il.next)
##      if (il.cur.info.degree > 1) COPYANDDEC(il.cur)
##      }
##  SWITCHCOORDS(a,b)
##  if (a.pa != noind)
##    {
##    if (a.pa.info.degree > 1) { DISTABOUT(a.pa) curtot += distab; }
##    if (a.ma.info.degree > 1) { DISTABOUT(a.ma) curtot += distab; }
##    for (il= a.kids; il != noindlist; il = il.next)
##      if (il.cur.info.degree > 1) { DISTABOUT(il.cur) curtot += distab; }
##    }
##  if (b.pa != noind)
##    {
##    if (b.pa.info.degree > 1) { DISTABOUT(b.pa) curtot += distab; }
##    if (b.ma.info.degree > 1) { DISTABOUT(b.ma) curtot += distab; }
##    for(il= b.kids;il != noindlist;il = il.next)
##      if (il.cur.info.degree > 1) { DISTABOUT(il.cur) curtot += distab; }
##    }
##
##  ''' Acceptance criterion '''
##  if (curtot < oldtot)
##    { lastdown = i; ndown++; }
##  else if (curtot == oldtot) 
##    {lastacross = i; nacross++; }
##  else if (curtot  <= oldtot - temp*log(randu01()))
##    { lastup = i; freezing = temp; nup++; }
##  else 
##    {
##    curtot = oldtot;
##    xt = a.info.xcoord;
##    yt = a.info.ycoord;
##    a.info.xcoord = b.info.xcoord;
##    a.info.ycoord = b.info.ycoord;
##    b.info.xcoord = xt;
##    b.info.ycoord = yt;
##    for (j = indpt + 1; j <= newpt; j++)
##      {
##      indstack[j].pa.info.xcoord = indstack[j].info.xcoord;
##      indstack[j].pa.info.ycoord = indstack[j].info.ycoord;
##      }
##    }
##  }
##
##''' End of main loop '''
##
##''' Output results '''
##fprintf(stderr,'\n Total squared x distances after searching= %12.5f',curtot);
##fprintf(stderr,'\n Finishing temperature                    = %12.5f',temp);
##fprintf(stderr,'\n Freezing temperature                     = %12.5f',freezing);
##fprintf(stderr,'\n Number of tries                          = %8d',ntries);
##fprintf(stderr,'\n Number of downhill steps                 = %8d',ndown);
##fprintf(stderr,'\n Last downhill step                       = %8d',lastdown);
##fprintf(stderr,'\n Number of uphill steps                   = %8d',nup);
##fprintf(stderr,'\n Last uphill step                         = %8d',lastup);
##fprintf(stderr,'\n Number of lateral steps                  = %8d',nacross);
##fprintf(stderr,'\n Last lateral step                        = %8d\n',lastacross);
##
##''' Set coords of all individuals '''
##for (i=1, ndiscon=1; i<=indpt; i++)
##  if (indstack[i].info.degree == 0) ndiscon++;
##  else if (indstack[i].info.degree < 2) setindcoords(&indstack[i],y_unit);
##for (i=1, j=0; i<=indpt; i++)
##  if (indstack[i].info.degree == 0) 
##    {
##    j++;
##    indstack[i].info.xcoord = 1.0*j/((float) ndiscon);
##    indstack[i].info.ycoord = .03;
##    }
##    
##for (i=1; i<=marpt; i++)
##  if (marstack[i].pa.par == nomar || marstack[i].ma.par == nomar)
##    if (marstack[i].pa.info.xcoord == marstack[i].ma.info.xcoord)
##      { 
##      if (marstack[i].pa.par == nomar) marstack[i].pa.info.xcoord -= spacer;
##      if (marstack[i].ma.par == nomar) marstack[i].ma.info.xcoord += spacer;
##      }
##}

def sort_index(a):
    '''Sort a list a by ascending value; return the corresponding list of indices a.'''
    return sorted(range(len(a)), key=a.__getitem__) 

'''
============================================================
Postscript pedigree drawing tools.

Corresponds to pedfiddler files:
postscript.h
============================================================
'''
########################################################################
# Constants
########################################################################
SQRT_2              = 0.707107
SQRT_PI             = 1.772454
UNITSPERPOINT       = 10
UNITSPERINCH        = 72 * UNITSPERPOINT
TITLESEP            = 1.0   # fraction of title point size
TICLEN              = 0.01   # fraction of plot size
DEFAULTINSIDE       = 0.04   # fraction of plot size -- room between box and plot
LABELSEP            = 0.40   # fraction of label point size
TWICE_MIN_DOT_SIZE  = 30.0   # multiple of line width

MEDIA = {
         'Letter'     : (  612,  792),
         'A4'         : (  595,  842),
         'Legal'      : (  612, 1008),
         'Statement'  : (  396,  612),
         'Tabloid'    : (  792, 1224),
         'Ledger'     : (  612,  792),
         'Folio'      : (  612,  936),
         'Quarto'     : (  610,  792),
         '10x14'      : (  720, 1008)
         }

########################################################################
# Public Methods
########################################################################
def draw_pedigree_to_file(ped, ped_info, plot_params, out_file):
    '''The main method that draws a pedigree in EPS format into the file name or file stream
    out_file using the plot parameters in plot_params.'''
    if isinstance(out_file, str):
        with open(out_file, 'wb') as out:
            plotter = __EpsPlotter(ped, ped_info, plot_params, out)
            plotter.draw_pedigree()
    else:
        plotter = __EpsPlotter(ped, ped_info, plot_params, out_file)
        plotter.draw_pedigree()

#========================================================================================
class PlotParams(object):
    def __init__(self):
        self.creator        = 'pedc2ps'
        self.username       = 'oren'
        self.title          = 'Pedigree'
        self.frame          = True
        self.ids            = True
        self.icons          = True
        self.labels         = True 
        self.rotate         = False 
        self.media          = None
        self.height         = 11.0
        self.width          = 8.5 
        self.linewidth      = 0.5   # in points
        self.margins        = 0.5
        self.titlefontsize  = 18.   # in points
        self.labelfontsize  = 9.    # in points
        self.iconsize       = 9.    # in points
        self.elbow_room     = 0.    # in points
        self.extra_margin   = 0.

########################################################################
# Private Methods
########################################################################
class __EpsPlotter(object):
    '''Plots a pedigree into an output stream .'''
    def __init__(self, ped, ped_info, plot_params, out):
        self._ped       = ped
        self._p         = plot_params
        self._out       = out
        self._ped_info  = ped_info.copy()
        self._box_coord = None
        self._box_node  = None

    ########################################################################
    # Public Methods
    ########################################################################
    def draw_pedigree(self):
        '''Draws the pedigree as a marriage node graph to the standard output file using the
        plot parameters in p and the node positions in ped_info.
           If ( p.frame ),   a frame is drawn.
           If ( p.icons ),   individuals are marked with icons. 
                             Else, individuals are not marked with icons.
           If ( p.title[0] ) p.title is the plot title, else no title appears 
        
           p.width is the width of the page in inches.
           p.height is the height of the page in inches 
             (including the title and any legend).
             The plot is centered in a region PAGEWIDTH by PAGEHEIGHT inches.
        
           p.titlefontsize is the type size for the title in points ( 1pt = 1/72 inch ).
           p.labelfontsize is the type size for all other text on the plot in points.
           p.iconsize is the size in points of the side the square icon (male).
           p.linewidth is the thickness of all lines
              (frame, graph links, icon outlines ) in points.
           p.elbow_room is the minimum separation between tokens (approx).
           p.extra_margin is the extra space in points between the plotted 
              individuals and marriages and the sides of the box -- default 0.0'''
        (p, ped, ped_info) = (self._p, self._ped, self._ped_info)
#        d = p if isinstance(p, dict) else p.__dict__
#        print sorted(d.iteritems(), key=operator.itemgetter(0))
        
        # Initialize plot
        self._box = self.__bounding_box(p)
        box = self._box
        self.__draw_header(p, box)
        self.__draw_title(p, box)
        (self._box_coord, self._box_node) = self.__create_margin(p, box, ped_info)
        
        #if p.elbow_room > 0:
        #    self.__rescale_x((UNITSPERPOINT*(p.iconsize+p.elbow_room))/((float)(box[1]-box[0])))
        #    ped.balance_marriages()
        #    (self._box_coord, self._box_node) = self.__create_margin(ped)
        if p.frame:
            self.__draw_frame_plot(box)
        
        # Draw data
        self.__draw_edges(ped, ped_info.data)
        self.__draw_nodes(p, ped, ped_info.data)
        
    ########################################################################
    # Private Methods
    ########################################################################            
    def __draw_nodes(self, p, ped, info_data):
        '''Plot person nodes.'''
        if p.icons:
            for name in ped.person_nodes:
                (node, info) = (ped.data[name], info_data[name])
                coord = self.__rel_coord(info)
                
                self.__rgb_icon(info.rgb)
                self.draw_node(node.sex, coord)
          
                if p.ids and info.name:
                    if (node.sex == SEX.UNKNOWN):
                        delta = ( LABELSEP * p.labelfontsize +  SQRT_2 * p.iconsize ) * UNITSPERPOINT
                    else:
                        delta = ( LABELSEP * p.labelfontsize + p.iconsize / SQRT_PI ) * UNITSPERPOINT
                    self.__draw_label((coord[0], coord[1] + delta), info.name)
          
                if p.labels and node.label:
                    if (node.sex == SEX.UNKNOWN):
                        delta = ( (0.75+LABELSEP)* p.labelfontsize + SQRT_2 * p.iconsize ) * UNITSPERPOINT
                    else:
                        delta = ( (0.75+LABELSEP)* p.labelfontsize + p.iconsize / SQRT_PI ) * UNITSPERPOINT
                    self.__draw_label((coord[0], coord[1] - delta), node.label)
        self.__ps_footer()

    def draw_node(self, sex, coord):
        if (sex == SEX.FEMALE):
            self.__draw_circle(coord)
        elif (sex == SEX.MALE):
            self.__draw_square(coord)
        elif (sex == SEX.UNKNOWN):
            self.__draw_diamond(coord)
        
    def __draw_edges(self, ped, info_data):
        '''Plot edges and marriage nodes.'''
        g = ped.undirected_graph
        for node in ped.nodes_with_type(NODE_TYPE.MARRIAGE):
            info = info_data[node]
            # Plot all family edges
            for nbhr in g.neighbors(node):
                self.__draw_edge(info, info_data[nbhr])
            # Plot marriage node
            self.__draw_dot_mar(info)
            self.write('Sk\n')
        
    def __draw_title(self, p, box):
        # TODO: add setprecision(4)
        if p.title:
            self.write('gsave\n'
                      + 'currentfont ' + repr(int(p.titlefontsize / p.labelfontsize)) 
                      + ' scalefont setfont\n'
                      + '(' + p.title + ') dup stringwidth pop 2.0 div '
                      + repr(box.llx + ( box.urx - box.llx ) / 2) + ' sub neg '
                      + repr(box.buy + box.titleroom / 2) + ' moveto show\n'
                      + 'grestore' + '\n')
                  
    def __bounding_box(self, p):
        box = _BoundingBox()
        box.titleroom = iround( UNITSPERPOINT * p.titlefontsize * ( 1.0 + TITLESEP ) ) if p.title else 0 
        margin = iround( UNITSPERINCH * p.margins )
        box.llx = margin
        box.lly = margin
        box.bly = box.lly
        box.urx = iround( UNITSPERINCH * (p.height if p.rotate else p.width)) - margin
        box.ury = iround( UNITSPERINCH * (p.width if p.rotate else p.height)) - margin
        box.buy = box.ury - box.titleroom
        return box
    
    def __create_margin(self, p, box, ped_info):
        maxxcoord = max(ped_info.coord(0))
        minxcoord = min(ped_info.coord(0)) 
        maxycoord = max(ped_info.coord(1))
        minycoord = min(ped_info.coord(1)) 
    
        tic_length = iround( TICLEN * min( (box.urx - box.llx, box.buy - box.bly) ) )
        '''onespace is the space from either side to any individual or marriage '''
        onespace = iround(p.iconsize/SQRT_PI * UNITSPERPOINT / 2 + 
                       DEFAULTINSIDE * min( (box.urx - box.llx, box.buy - box.bly) ))
        '''onespace + otherspace it the space form the top or bottom to 
             any individual or marriage '''
        otherspace = iround((1.50 + LABELSEP) * p.labelfontsize * UNITSPERPOINT)
        '''plx = llx + tic_length + onespace, prx = urx - onespace'''
        plx = iround(box.llx + tic_length + onespace + p.extra_margin * UNITSPERPOINT)
        prx = iround(box.urx - onespace - p.extra_margin * UNITSPERPOINT)
        ply = iround(box.bly + onespace + otherspace)
        puy = iround(box.buy - onespace - otherspace)
        
        return ((minxcoord, maxxcoord, minycoord, maxycoord), (plx, prx, ply, puy))

    def __draw_header(self, p, box):
        csze = int(UNITSPERPOINT * p.iconsize / SQRT_PI)
        ssze = int(UNITSPERPOINT * p.iconsize)
        asze = int(UNITSPERPOINT * p.iconsize * SQRT_2)
        dsze = int(max( (TWICE_MIN_DOT_SIZE * p.linewidth, csze * p.linewidth * 0.75) ))
                  
        cex  = p.labelfontsize / 14.0
        lwd  = int( 40.0 * cex )
    
        medN = p.media
        (medX, medY) = MEDIA[medN] if medN else (None, None)
    
        self.write('%!PS-Adobe-3.0')
        if p.media:
            self.write('\n')
        else:
            self.write(' EPSF-3.0\n'
                       + '%%DocumentNeededResources: fonts Helvetica\n'
                       + '%%DocumentFonts: Helvetica\n')
        if not p.title:
            self.write('%%Title: (none)\n')
        else:
            self.write('%%Title: ' +  p.title + '\n'
                       + '%%Creator: ' + p.creator + '\n'
                       + '%%CreationDate: ' +  time.ctime(time.time()) + '\n'
                       + '%%For: ' + p.username + '\n'
                       + '%%DocumentData: Clean7Bit' + '\n')
        if p.media:
            self.write('%%DocumentMedia: ' + medN + ' ' + medX + ' ' + medY 
                       + ' 75 ( ) ( )\n'
                       +  '%%Pages: 1\n')
        bbf = int(3+p.linewidth) if p.frame else 0 
        if ((int(box.llx / 10.0)) < bbf or (iround(box.lly / 10.0)) < bbf):
            bbf = 0
        self.write('%%BoundingBox: '
                   + repr(int(box.llx / 10.0) - bbf) + ' ' 
                   + repr(int(box.lly / 10.0) - bbf) + ' '
                   + repr(int((box.ury if p.rotate else box.urx) / 10.0) + bbf) + ' '
                   + repr(int((box.urx if p.rotate else box.ury) / 10.0) + bbf) + '\n'
                   + '%%EndComments\n'
                   + '%%BeginProlog: procset Pedfiddler_proc\n'
                   + '/Pedfiddler_proc 100 dict def\n'
                   + 'Pedfiddler_proc begin\n'
                   + '/bd {bind def} def\n'
                   + '/M {moveto} bd\n'
                   + '/L {lineto} bd\n'
                   + '/R {rlineto} bd\n'
                   + '/C {closepath} bd\n'
                   + '/N {newpath} def\n'
                   + '/Sk {stroke} def\n'  
                   + '/GS {gsave} def\n'
                   + '/GR {grestore} def\n'
                   + '/SGF {setgray fill} bd\n'
                   + '/tex {\n'
                   + '  GS\n'
                   + '  M 0 index\n'
                   + '  dup stringwidth pop neg 2 div 0 rmoveto\n'
                   + '  GS\n  false charpath\n  1 setgray ' + repr(lwd) + ' setlinewidth\n'
                   + '  1 setlinejoin Sk GR\n'
                   + '  show GR\n'
                   + '  } bd\n'
                   + '/edge {\n  M L\n} bd\n'
                   + '/dot {\n'
                   + '  GS N ' + repr(dsze) +' 0 360 arc 0 SGF GR\n'
                   + '  } bd\n'
                   + '/cir {\n'
                   + '  N ' + repr(csze) + ' 0 360 arc\n'
                   + '  GS setrgbcolor fill GR Sk\n'
                   + '  } bd\n'
                   + '/squ {\n'
                   + '  N M\n'
                   + '  0 ' + repr(ssze) + ' R ' + repr(ssze) + ' 0 R\n'
                   + '  0 ' + repr(-ssze) + ' R ' + repr(-ssze) + ' 0 R\n'
                   + '  C\n'
                   + '  GS setrgbcolor fill GR Sk\n'
                   + '  } bd\n'
                   + '/dia {\n'
                   + '  N M\n'
                   + '  ' + repr(-asze) + ' ' + repr(asze) + ' R ' + repr(asze)
                   + ' ' + repr(asze) +' R\n'
                   + '  ' + repr(asze) + ' ' + repr(-asze) + ' R ' + repr(-asze) 
                   + ' ' + repr(-asze) +' R\n'
                   + '  C\n'
                   + '  GS setrgbcolor fill GR Sk\n'
                   + '  } bd\n'
                   + 'end\n'
                   + '%%EndProlog\n'
                   + '%%BeginSetup' + '\n')
        if p.media:
            self.write('%%BeginFeature: *PageSize: ' + repr(medN) + '\n'
                       + '  + /PageSize [' + repr(medX) + ' ' + repr(medY) 
                       + '] /ImagingBBox null >> setpagedevice\n'
                       + '%%EndFeature' + '\n')
        self.write('Pedfiddler_proc begin\n'
                   + 'gsave\n'
                   + '%%EndSetup\n'
                   + '%%Page: 1 1\n'
                   + repr(iround(UNITSPERPOINT * p.linewidth)) + ' setlinewidth\n'
                   + '0.1 0.1 scale\n'
                   + '1 setlinecap 1 setlinejoin\n'
                   + '0 setgray\n'
                   + '%%IncludeResource: font Helvetica\n'
                   + '/Helvetica findfont ' + repr(int(140*cex)) 
                    + ' scalefont setfont' + '\n')
        if p.rotate:
            self.write(repr(iround(p.width * UNITSPERINCH)) + ' 0 translate 90 rotate\n')
        
    def __rel_coord_of_index(self, x, i):
        box = self._box_node
        (l, r) = (2*i, 2*i+1)
        (min_coord, max_coord) = self._box_coord[l:r+1]
        return (box[l]+(box[r]-box[l])*(x-min_coord)/(max_coord-min_coord)) if max_coord-min_coord > 0 \
            else 0.5*(box[l]+box[r])
    
    def __rel_coord(self, m):
        return (self.__rel_coord_of_index(m.coord[0], 0), self.__rel_coord_of_index(m.coord[1], 1))
    
    def __draw_edge(self, m, i):
        self.write('%d %d %d %d edge\n' % (self.__rel_coord(m) + self.__rel_coord(i)))

    def __draw_dot_mar(self, m):
        self.write('%d %d dot\n' % self.__rel_coord(m))

    def __draw_frame_plot(self, box):
        '''Make box with nx tick marks on x-axis and ny tick marks on y-axis.'''
        self.write('N\n'
                   + repr(box.llx) + ' ' + repr(box.bly) + ' M\n'
                   + repr(box.llx) + ' ' + repr(box.buy) + ' L\n'
                   + repr(box.urx) + ' ' + repr(box.buy) + ' L\n'
                   + repr(box.urx) + ' ' + repr(box.bly) + ' L\n'
                   + repr(box.llx) + ' ' + repr(box.bly) + ' L\n'
                   + 'C Sk\n')
                   
    def __rgb_icon(self, rgb):
        self.write('%.4f %.4f %.4f ' % (rgb[0]/255.0, rgb[1]/255.0, rgb[2]/255.0))
    
    def __draw_square(self, x):
        r = UNITSPERPOINT * self._p.iconsize / 2.0
        self.write('%d %d squ\n' % (int(x[0]-r), int(x[1]-r)))
    
    def __draw_diamond(self, x):
        r = UNITSPERPOINT * self._p.iconsize * SQRT_2
        self.write('%d %d dia\n' % (int(x[0]-r), int(x[1]-r)))
    
    def __draw_circle(self, x):
        self.write('%d %d cir\n' % x)
    
    def __draw_label(self, x, s):
        self.write('(%s) %d %d tex\n' % ((s,) + x))
    
    def __ps_footer(self):
        self.write('showpage\n'
                   + '%%Trailer\n'
                   + 'grestore\n'
                   + 'end\n'
                   + '%%EOF' + '\n')

    def write(self, s):
        self._out.write(s)

#========================================================================================
class _BoundingBox(object):
    '''BoundingBox is 
       llx (lower left x coord)
       lly (lower left y coord)
       urx (upper left x coord)
       ury (upper left y coord)
    
       bly = lly : y coord of lower edge of box
       buy       : y coord of upper edge of box
    
       box_node:
       plx : minimum x coord of individual or marriage
       prx : maximum x coord of individual or marriage
       ply : minimum y coord of individual or marriage
       puy : maximum y coord of individual or marriage
    '''
    def __init__(self):
        self.llx = None
        self.lly = None
        self.urx = None
        self.ury = None

        self.bly = None
        self.buy = None

        self.plx = None
        self.prx = None
        self.ply = None
        self.pry = None

        self.titleroom = None

#========================================================================================
    
#    def __rescale_x(self, ped_info, delta):
#        '''Rescale margins to accomodate data.'''
#        maxxcoord = max(ped_info.coord(0))
#        minxcoord = min(ped_info.coord(0)) 
#        xdelta = delta * ( maxxcoord - minxcoord ) if (maxxcoord - minxcoord > 0) else 1
#        
#        indpt = len(indstack)
#        ngen = [0] * indpt
#    
#        k = 0
#        ycoord = ped_info.coord(1)
#        for (i, y) in enumerate(ycoord):
#            for j in xrange(0, i):
#                if (y[j] == y):
#                    k += 1
#                    break
#            
#        gen = indpt - k
#        idx = [0] * gen
#        idxbuf = [0] * indpt
#        lastycoord = - 10000.0
#        left = indpt
#        j = 0
#    
#        while left > 0:
#            idx[j] = idxbuf + indpt - left
#            ycoord = 10000.0
#            for i in xrange(1, indpt+1):
#                if (info.ycoord > lastycoord and info.ycoord < ycoord ):
#                    ycoord = info.ycoord
#            ngen[j] = 0
#            for i in xrange(1, indpt+1):
#                if(info.ycoord == ycoord):
#                    idx[j][ngen[j]] = i
#                    ngen[j] += 1
#                    left -= 1
#            j += 1
#            lastycoord = ycoord
#      
#        '''Sort generations'''
#        for i in xrange(0, gen):
#            for j in xrange(0, ngen[i]-1):
#                keyindex = idx[i][j]
#                keycoord = indstack[keyindex].info.xcoord
#                for k in xrange(j+1, ngen[i]):
#                    if (keycoord > indstack[idx[i][k]].info.xcoord):
#                        keyindex = idx[i][k]
#                        keycoord = indstack[keyindex].info.xcoord
#                        idx[i][k] = idx[i][j]
#                        idx[i][j] = keyindex
#            idx[i][j] = keyindex
#    
#        npools = [0] * gen
#        pools = [0] * gen
#        poolsbuf = [0] * indpt
#    
#        offset = 0
#        for i in xrange(0, gen):
#            pools[i] = poolsbuf + offset
#            offset += ngen[i]
#    
#        poolsum = [0.0] * gen
#        poolsumbuf = [0.0] * indpt
#    
#        offset = 0
#        for i in xrange(0, gen):
#            poolsum[i] = poolsumbuf + offset
#            offset += ngen[i]
#    
#        for i in xrange(0, gen):
#            npools[i] = 0
#            for j in xrange(0, ngen[i]):
#                npools[i] += 1
#                pools[i][npools[i]-1] = 1
#                poolsum[i][npools[i]-1] = indstack[idx[i][j]].info.xcoord
#                while(npools[i] > 1 and
#                      ((poolsum[i][npools[i]-2] / pools[i][npools[i]-2]) +
#                       (pools[i][npools[i]-2] + pools[i][npools[i]-1] ) * xdelta * 0.5
#                       > ( poolsum[i][npools[i]-1] / pools[i][npools[i]-1]))):
#                    npools[i] -= 1
#                    pools[i][npools[i]-1] += pools[i][npools[i]]
#                    poolsum[i][npools[i]-1] += poolsum[i][npools[i]]
#    
#        for i in xrange(0, gen):
#            j = 0
#            for jj in xrange(0, npools[i]):
#                xstart = poolsum[i][jj] / pools[i][jj] - ( pools[i][jj] - 1 ) *  xdelta * 0.5
#                if ( pools[i][jj] > 1 ):
#                    for k in xrange(0, pools[i][jj]):
#                        indstack[idx[i][j]].info.xcoord = xstart + k * xdelta
#                        j += 1
#                else:
#                    j += 1
#
#    def __sort_generations(self, y):
#        [y_sorted, index] = sorted((y, range(len(y))), key = lambda x: x[0])
#        gen_start_index = __group_index(index)
#
#def __group_index(a):
#    '''In a collection a, return the indices of the start of each group of identical
#    consecutive elements.'''
#    (i_prev, x_prev) = (-1, None)
#    for (i, x) in enumerate(a):
#        if i_prev < 0 or x_prev != x:
#            (i_prev, x_prev) = (i, x)
#            yield i
            
def iround(num):
    '''Round to the nearest integer.'''
    return int(num + .5) if num > 0 else int(num - .5)

########################################################################
# Public Methods
########################################################################
def draw_pedigree(in_file, coord_params, plot_params, colors, out_file):
    '''Main pedigree drawing function.'''
    ped = read_pedigree(in_file)
    ped_info = PedigreeInfo(ped)
    compute_coords(ped_info, coord_params)
    ped_info.set_rgb(colors)
    draw_pedigree_to_file(ped, ped_info, plot_params, out_file)
    return ped_info

########################################################################
# Argument Parsing
########################################################################
def __parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    prog = os.path.basename(sys.argv[0])
    usage = 'Usage: %s [options] <in_file> <out_file>\n\n' \
        'Draw a pedigree. This is a python translation of the pedfiddler algorithm\n' \
        'by Loredo-Osti and Morgan (http://www.stat.washington.edu/thompson/Genepi/Pedfiddler.shtml)\n\n' \
        'Input: a text file whose row format is\n' \
        '<id> <father_id> <mother_id> <sex> <affection_status> <label>\n\n' \
        'ID coding: missing=0, regular ID=any positive number\n' \
        'Sex coding: unknown=0, male=1, female=2\n' \
        'Affection status coding:  unknown=0, normal=1, affected=2\n\n' \
        'Output: encapsulated postscript (EPS) file\n\n' \
        'Type ''%s -h'' to display full help.' % (prog, prog)

    parser = OptionParser(usage=usage, option_class=MyOption)
    parser.add_option('-t', '--title'        , default='Pedigree', dest='title', help='Figure title')
    parser.add_option('-f', '--frame'        , action='store_true', dest='frame', default=True, help='Draw frame [boolean]')
    parser.add_option('-i', '--ids'          , action='store_true', dest='ids', default=True, help='Draw node IDs [boolean]')
    parser.add_option('-c', '--icons'        , action='store_true', dest='icons', default=True, help='Draw node icons [boolean]')
    parser.add_option('-l', '--labels'       , action='store_true', dest='labels', default=True, help='Draw node labels [boolean]')
    parser.add_option('-r', '--rotate'       , action='store_true', dest='rotate', default=False, help='Rotate figure [boolean]')
    parser.add_option('-p', '--media'        , dest='media', default=None, help='EPS paper type (%s)' % '/'.join(MEDIA.iterkeys()))
    parser.add_option('-H', '--height'       , type="float", default=11.0, help='Figure height [inches]')
    parser.add_option('-W', '--width'        , type="float", default=8.5, help='Figure weight [inches]')
    parser.add_option('-w', '--linewidth'    , type="float", default=1.0, help='Edge line width [points]')
    parser.add_option('-m', '--margins'      , type="float", default=0.5, help='Edge line width [points]')
    parser.add_option('-T', '--titlefontsize', type="float", default=18.0, help='Title font size [points]')
    parser.add_option('-L', '--labelfontsize', type="float", default=9.0, help='Label font size [points]')
    parser.add_option('-I', '--iconsize'     , type="float", default=14.0, help='Icon font size [points]')
    parser.add_option('-v', '--debug'        , action='store_true', dest='debug', default=False, help='Print debugging information')
    parser.add_option('-U', '--unknown'      , type='rgb', dest = 'color_unknown', default=DEFAULT_COLORS[STATUS.UNKNOWN], help='Affected individual color [e.g., (255,255,255)]')
    parser.add_option('-N', '--normal'       , type='rgb', dest = 'color_normal', default=DEFAULT_COLORS[STATUS.NORMAL], help='Affected individual color [e.g., (255,255,255)]')
    parser.add_option('-A', '--affected'     , type='rgb', dest = 'color_affected', default=DEFAULT_COLORS[STATUS.AFFECTED], help='Affected individual color [e.g., (255,255,255)]')
    (options, args) = parser.parse_args(sys.argv[1:])
    if len(args) != 2:
        print usage
        sys.exit(1)
    (options.input, options.output) = args
    (valid, options) = __validate_options(options)
    if not valid:
        print usage
        sys.exit(1)
    return options
    
def check_rgb(option, opt, value):
    '''Parse an RGB tuple string like ''(255,255,0)'' to the tuple (255,255,0).'''
    try:
        return reduce(tuple.__add__, int(value[1:-1].split(',')))
    except ValueError:
        raise OptionValueError(
            "option %s: invalid RGB tuple value: %r" % (opt, value))    

class MyOption(Option):
    TYPES = Option.TYPES + ("rgb",)
    TYPE_CHECKER = deepcopy(Option.TYPE_CHECKER)
    TYPE_CHECKER["rgb"] = check_rgb
    
def __validate_options(options):
    '''Validate input options. Set default/overrides.'''
    # Set defaults/overrides
    if not os.path.isfile(options.input):
        return (False, options)

    # Coordinate computation options
    options.algorithm = 'default'
    options.balance_marriages = True
    #options.num_iterations = 1000
    #options.initial_temp = 15.0
    #options.cooling_factor = 0.9    
    options.colors = {STATUS.UNKNOWN : options.color_unknown,
                      STATUS.NORMAL  : options.color_normal,
                      STATUS.AFFECTED: options.color_affected}
    # Plot options
    options.creator = 'famplot'
    options.username = 'username'
    options.elbow_room = 0.0
    options.extra_margin = 0.0

    return (True, options)

########################################################################
# Main Program
########################################################################
if __name__ == '__main__':
    '''Main program interface - accepts CLI arguments.'''
    options = __parse_command_line_args()
    if options.debug:
        print 'Input options:', options
    draw_pedigree(options.input, options, options, options.colors, options.output)
