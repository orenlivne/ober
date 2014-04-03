'''
============================================================
A pedigree of haploids (e.g., humans). This class is not
thread-safe. See [1] for information on properly handling the
lazily-initialized depth field. It contains the basic data
structures and functions related to the pedigree; see PetTools
for manipulative tools.

This object is meant to be immutable.

Created on May 31, 2012
@author: Oren Livne <livne@uchicago.edu>

@see [1] http://michelanders.blogspot.com/2010/12/python-thread-safe-cache-class.html
============================================================
'''
import numpy as np, networkx as nx, util, impute as im
from networkx.algorithms.dag import is_directed_acyclic_graph

class Pedigree(object):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    # ID to start generating IDs when recoding IDs
    START_ID = 0  # Since we use it also as a python array index
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, graph, sample_id=None, sex=None, phenotype=None, node_type=None,
                 sample_index=None, num_genotyped=0):
        '''Construct a pedigree from a Person DAG (edge goes from parent to child).'''
        # Validate input
        if not is_directed_acyclic_graph(graph): raise ValueError('Pedigree graph must be a directed acyclic graph (DAG)')
        for node in graph.nodes_iter():
            degree = graph.in_degree(node)
            if degree > 2: raise ValueError('In-degree must be at most 2, but degree[%d] = %d' % (node, degree,))
 
        # Set defaults
        if sample_id is None:
            sample_id = graph.nodes()
            sex = np.tile(im.constants.INDETERMINATE, (graph.number_of_nodes(), 1))
            node_type = np.tile(im.constants.INDETERMINATE, (graph.number_of_nodes(), 1))
            phenotype = np.tile(im.Person.Person.PHENOTYPE_MISSING, (graph.number_of_nodes(), 1))

        # Set input fields
        self._graph = graph
        self._graph_undirected = graph.to_undirected()
        self._sample_id = sample_id
        self._sample_index = sample_index if sample_index is not None else (np.arange(0, len(sample_id)) if sample_id is not None else None)
        self._num_genotyped = num_genotyped
        
        # Set derived fields
        # Reverse map sample-id -> data index if ids have been recoded
        self._person = dict((node, im.Person.Person(node, sex[index], phenotype[index], node_type[index],
                                                    self.father(node, True), self.mother(node, True)))
                            for (index, node) in enumerate(graph.nodes_iter()))
        self._node_of = dict(zip(sample_id   , range(0, self.n))) if sample_id    is not None else None
        self._index_of = dict(zip(sample_index, range(0, self.n))) if sample_index is not None else None

        # Lazily-initialized cached properties - node ordering        
        self._depth = None
        self._preorder = None
        self._postorder = None

        # Lazily-initialized cached properties        
        self._duos = {}
        self._trios = {}
        self._families = {}
        self._snp_range = None


    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __key(self):
        '''Uniquely-identifying key of this object.'''
        pass
    
    def __eq__(self, y):
        '''Equality of objects.'''
        return self.__key() == y.__key()

    def __ne__(self, y):
        '''Inequality of objects.'''
        return self.__key() != y.__key()

    def __hash__(self):
        '''Object hash code.'''
        return hash(self.__key())

    def __repr__(self):
        return 'Pedigree[n=%d, samples=%d, genotyped=%d]' % (self.n, self.n, self.num_genotyped) 
    
    ####################################################################################
    # IDs and finding objects by IDs
    ####################################################################################
    def sample_ids(self, node):
        '''Return the sample IDs of the node indices in the collection ''node''.'''
        return np.array([self._sample_id[x] for x in node]) 
    
    def genotyped_sample_id(self):
        return np.array([self._sample_id[k] for k, _ in enumerate(self.sample_index) if self.is_genotyped(k)])

    def genotyped_sample_index(self):
        return np.array([x for k, x in enumerate(self.sample_index) if self.is_genotyped(k)])
    
    def person_with_ids(self, i):
        '''a Person object corresponding to node i with ID, father ID and mother ID instead of indices.'''
        person = self.person[i].copy()
        (person.sample_id, person.father, person.mother) = self.trio_sample_id(i)
        return person
    
    def is_genotyped(self, x):
        '''Does the sample index x correspond to a genotyped sample or not.''' 
        return (x < self._num_genotyped)

    ####################################################################################
    # Parents and children
    ####################################################################################
    def parents(self, i):
        '''Return a dictionary of parent type (PATERNAL/MATERNAL) to parent node ID of node i.'''
        parent_nodes = self._graph.predecessors(i)
        parent_types = [self._graph[p][i]['type'] for p in parent_nodes]
        return dict(zip(parent_types, parent_nodes)) 

    def parent(self, i, parent_type, nullsafe=False):
        '''Return the parent of node i of type parent_type (PATERNAL/MATERNAL).'''
        parents = self.parents(i)
        return (parents[parent_type] if parents.has_key(parent_type) else im.constants.MISSING) if nullsafe else parents[parent_type] 

    def father(self, i, nullsafe=False):
        '''Return the father of node i.'''
        return self.parent(i, im.constants.PATERNAL, nullsafe=nullsafe)

    def mother(self, i, nullsafe=False):
        '''Return the mother of node i.'''
        return self.parent(i, im.constants.MATERNAL, nullsafe=nullsafe)

    def trio(self, i):
        '''Return the trio of i and its parent node IDs. Note: will return less than 3 elements
        in the tuple if i's parents are not in the pedigree.'''
        return (i, self._graph.predecessors(i))
        
    def neighbors(self, node, depth=1, genotyped=False, add_parents=False):
        '''Return all neighbors of a node in the entire pedigree (if genotyped = False) or
        all genotyped neighbors (if genotyped = True) of depth <= depth.
        If add_parents=False, adds all parents of all neighbors to complete the pedigree.'''
        nbhrs = [x for x in nx.single_source_shortest_path(self._graph_undirected, node, cutoff=depth).iterkeys() 
                 if (self.is_genotyped(x) if genotyped else True)]
        if add_parents: nbhrs = list(set(nbhrs + [parent for x in nbhrs for parent in self.parents(x).itervalues()])) 
        return nbhrs

    def neighbors_genotyped_selective(self, node, depth=1):
        '''Return all neighbors of a node in the entire pedigree (if genotyped = False) or
        all genotyped neighbors (if genotyped = True) of depth <= depth. Selects only neighbors
        that have genotyped neighbors on the path between the neighbor and node.'''
        return list(reduce(set.union,
                           (v for (k, v) in 
                            nx.single_source_shortest_path(self._graph_undirected, node, cutoff=depth).iteritems() 
                            if self.is_genotyped(k)), set([])))
        
    def trio_sample_id(self, i):
        '''Return the sample ID tuple (ID, father ID, mother ID) of node i.'''
        ids, parents = self._sample_id, self.parents(i)
        return ids[i], \
            ids[parents[im.constants.PATERNAL]] if parents.has_key(im.constants.PATERNAL) else im.constants.MISSING, \
            ids[parents[im.constants.MATERNAL]] if parents.has_key(im.constants.MATERNAL) else im.constants.MISSING

    ####################################################################################
    # Duos, Trios & Family Lists
    ####################################################################################
    def duos(self, parent_type, genotyped=True):
        '''Find and cache all genotyped duos in the data set:
        (child-parent[parent_type]) tuples. Cached.'''
        return self._duos[parent_type] if self._duos.has_key(parent_type) else \
            self._duos.setdefault(parent_type, np.array(list(self.__compute_duos(parent_type, genotyped))))
    
    def trios(self, genotyped=True):
        '''Find and cache all genotyped trios in the data set. Sorted by parents. Cached.'''
        return self._trios[genotyped] if self._trios.has_key(genotyped) else \
            self._trios.setdefault(genotyped, np.array(list(Pedigree.__sorted_trios(self.__compute_trios(genotyped)))))

    def kids_duos(self, kids, parent_type, genotyped=True):
        '''Find and cache all genotyped duo IDs in the data set:
        (child-parent[parent_type]) tuples where child is in kids. Not cached.'''
        return np.array(list(self.__compute_duos(parent_type, genotyped, kids=kids)))

    def kids_trios(self, kids, genotyped=True):
        '''Return the trios of the kids in the list ''kids''. Not cached.'''
        return np.array(list(Pedigree.__sorted_trios(self.__compute_trios(genotyped, kids=kids))))

    def families(self, genotyped=True, min_children=0, max_children=np.inf):
        '''Return an iterator over families in a list of trios in which both parents are in G and at
        least min_children of their children are in G. min_children=0 by default. Cached.'''
        key = (genotyped, min_children)
        return self._families[key] if self._families.has_key(key) else \
            self._families.setdefault(key, [family for family in Pedigree.__trios_to_families(self.trios(genotyped)).itervalues() if family.num_children >= min_children and family.num_children <= max_children])
    
    def families_union(self, genotyped=True, min_children=0, max_children=np.inf): 
        '''Return the union of all genotyped nuclear family members with at least min_children children
        in the trio list trio.'''
        try:
            return reduce(set.union, (family.member_set
                          for family in Pedigree.__trios_to_families(self.trios(genotyped)).itervalues() 
                          if family.num_children >= min_children and family.num_children <= max_children))
        except TypeError:
            # Found no such families => input to reduce() is empty
            return set([])
    
    ####################################################################################
    # Family Info and Selection
    ####################################################################################
    def find_family(self, father, mother, genotyped=True):
        '''Find a family by parents.'''
        f = [f for f in self.families(genotyped=genotyped) 
             if f.parents[im.constants.PATERNAL] == father and f.parents[im.constants.MATERNAL] == mother]
        return f[0] if f else None
        
    def find_families_by_either_parent(self, parent, genotyped=True, min_children=0):
        '''find_families_by_parent(p, parent_type, parent, min_children=0)
        Find families by parent type + ID.'''
        return [f for f in self.families(genotyped=genotyped, min_children=min_children) 
                if parent in f.parents]
        
    def find_families_by_parent(self, parent_type, parent, genotyped=True, min_children=0):
        '''find_families_by_parent(p, parent_type, parent, min_children=0)
        Find families by parent type + ID.'''
        return [f for f in self.families(genotyped=genotyped, min_children=min_children) 
                if f.parents[parent_type] == parent]
        
    def find_families_by_father(self, father_id, genotyped=True, min_children=0):
        '''Find families by father ID.'''
        return self.find_families_by_parent(im.constants.PATERNAL, father_id,
                                            genotyped=genotyped, min_children=min_children)
        
    def find_families_by_mother(self, mother_id, genotyped=True, min_children=0):
        '''Find families by mother ID.'''
        return self.find_families_by_parent(im.constants.MATERNAL, mother_id,
                                            genotyped=genotyped, min_children=min_children)

    def find_family_by_child(self, child_id, genotyped=True, min_children=0):
        '''Find the family that contain a child ID, or None, if not found.'''
        families = [f for f in self.families(genotyped=genotyped, min_children=min_children) if 
                    child_id in f.children]
        return families[0] if families else None

    def find_families_by_member(self, member_id, genotyped=True, min_children=0):
        '''Find all families that contain a member ID.'''
        return (f for f in self.families(genotyped=genotyped, min_children=min_children) if 
                member_id in f.member_set)

    def find_families_by_members(self, member_ids, genotyped=True, min_children=0):
        '''Find all families that contain a member ID.'''
        members = set(member_ids)
        return (f for f in self.families(genotyped=genotyped, min_children=min_children) if 
                members & f.member_set)

    ####################################################################################
    # Founders and Quasi-Founders
    ####################################################################################
    @property
    def quasi_founders(self):
        '''Return the array of all quasi-founder indices (genotyped samples whose parents are both
        not genotyped).''' 
        return np.where([all((y >= self.num_genotyped) for y in self.graph.predecessors(x)) for x in xrange(self.num_genotyped)])[0]

    @property
    def non_quasi_founders(self):
        '''Return the array of all quasi-founder indices (genotyped samples for which at least one
        parent is genotyped).''' 
        return np.setdiff1d(xrange(self.num_genotyped), self.quasi_founders)

    ####################################################################################
    # Misc
    ####################################################################################
    def sub_pedigree(self, samples):
        '''Return a sub-pedigree with the specified sample array only.'''
        if isinstance(samples, list): samples = np.array(samples)
        num_samples = len(samples)
        genotyped = np.where(samples < self.num_genotyped)[0]
        num_genotyped = len(genotyped) 
        samples = samples[np.concatenate((genotyped, np.where(samples >= self.num_genotyped)[0]))]
        
        # Create family pedigree graph
        # Original family IDs are self.sample_id[i] if ever needed
        nodes = range(0, num_samples) 
        recoding = dict(zip(samples, nodes))
        g = self.graph.subgraph(samples)
        sub_graph = nx.DiGraph()
        sub_graph.add_nodes_from(nodes)
        for parent_type in im.constants.ALLELES:
            sub_graph.add_edges_from(((recoding[n1], recoding[n2]) for n1, n2 in g.edges_iter() 
                                      if g.edge[n1][n2]['type'] == parent_type), type=parent_type)
        family_index = util.dict_invert(recoding)
        sample_index = [family_index[j] for j in nodes]
        sample_id = np.array([self.sample_id[family_index[j]] for j in nodes])
        sex = np.array([self.sex[family_index[j]] for j in nodes])
        phenotype = np.array([self.phenotype[family_index[j]] for j in nodes])
        node_type = np.array([self.node_type[family_index[j]] for j in nodes])
        
        return Pedigree(sub_graph, sample_id=sample_id, sex=sex, phenotype=phenotype, node_type=node_type,
                        sample_index=sample_index, num_genotyped=num_genotyped)
        
    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def graph(self):
        '''Pedigree A-directed Cyclic Graph (DAG) (edge goes from parent to child).'''  
        return self._graph
    
    @property
    def graph_undirected(self):
        '''Return the undirected graph counterpart.'''  
        return self._graph_undirected

    @property
    def person(self):
        '''a Person object list corresponding to graph nodes.'''  
        return self._person

    @property
    def n(self):
        '''Pedigree size (number of nodes).'''  
        return self._graph.number_of_nodes()
    
    @property
    def sample_id(self):
        '''# An array of sample_ids.'''  
        return self._sample_id

    @property
    def sample_index(self):
        '''# An array of sample_indices.'''  
        return self._sample_index
    
    @property
    def sample_id_genotyped(self):
        '''Return the sub-array of sample_ids of genotyped individuals.'''
        return self._sample_id[0:self._num_genotyped]

    @property
    def sex(self):
        '''Node sexes, as numpy array.'''  
        return np.array([self._person[i].sex for i in self._graph.nodes_iter()])

    @property
    def phenotype(self):
        '''Node phenotypes, as numpy array.'''  
        return np.array([self._person[i].phenotype for i in self._graph.nodes_iter()])

    @property
    def node_type(self):
        '''Node types, as numpy array.'''  
        return np.array([self._person[i].node_type for i in self._graph.nodes_iter()])

    @property
    def sex_dict(self):
        return dict((self._sample_id[i], self._person[i].sex) for i in self._graph.nodes_iter())

    @property
    def phenotype_dict(self):
        return dict((self._sample_id[i], self._person[i].phenotype) for i in self._graph.nodes_iter())
    
    @property
    def node_type_dict(self):
        return dict((self._sample_id[i], self._person[i].node_type) for i in self._graph.nodes_iter())
        return (self.graph.nodes(), set(self.graph.edges()), self.sample_id.tolist(),
                self.sex_dict, self.phenotype_dict, self.node_type_dict)
    
    def id_in_pedigree(self, node):
        return self._node_of.has_key(node)
        
    @property
    def num_genotyped(self):
        '''# gentyped samples in the sample_id map.'''  
        return self._num_genotyped

    @property
    def node_of(self):
        '''# Map sample-id -> node ID.'''  
        return self._node_of

    @property
    def index_of(self):
        '''# Map sample-id -> node index.'''  
        return self._index_of
    
    @property
    def is_founder(self):
        '''A flag array indicating whether a person is a founder or not..'''  
        return [not self._graph.predecessors(i) for i in self._graph.nodes_iter()] 
    
    @property
    def depth(self):
        '''Return the node depth ordering array. Once set, it is cached.'''
        return self._depth

    @depth.setter
    def depth(self, depth):
        '''Set a node depth ONLY if it is not set yet.'''
        if self._depth is None:
            self._depth = depth
        return self._depth

    @property
    def preorder(self):
        '''Create a pre-ordering of nodes on the problem's pedigree (parents before children) and cache.'''
        if self._preorder is None:
            # Use our generation numbers; networkx seems to produce wrong orders with dfs_preorder_nodes 
            depth = self.depth
            self._preorder = np.array(sorted(depth, key=depth.get))
        return self._preorder
    
    @property
    def postorder(self):
        '''Create a post-ordering of nodes on the pedigree (children before parents) and cache.'''
        if self._postorder is None:
            self._postorder = np.array(list(reversed(self.preorder)))
        return self._postorder
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __sample_index(self, genotyped):
        '''Return the sample indices of all samples or genotyped samples.'''
        # return self._sample_id if genotyped else xrange(0, self.n)
        return xrange(0, self.num_genotyped if genotyped else self.n)

    def __compute_duos(self, parent_type, genotyped=True, kids=None):
        '''Return a list of (x, parent[x]) where x=child ID and parent[x]=parent ID such that
        both the parent and child are genotyped.'''
        ids = self.__sample_index(genotyped)
        for i in (kids if kids is not None else ids):
            parents = self.parents(i)
            if parents.has_key(parent_type) and (util.is_member(ids, [parents[parent_type]]) if genotyped else True):
                yield (i, self.parents(i)[parent_type])

    def __compute_trios(self, genotyped=True, kids=None):
        '''Emit all (if kids=None) or selected genotyped trios for a pedigree p and genotype set g. Trios are unsorted.
        This assumes that parents are sorted in a consistent order (say, (father, mother)) in
        the lists returned from Pedigree.parents(). If kids is not None, returns only the trios of the kids in the list ids.'''
        ids = self.__sample_index(genotyped)
        for i in (kids if kids is not None else ids):
            parents = self.parents(i)
            if len(parents) == 2 and (util.is_member(ids, parents.itervalues()) if genotyped else True):
                yield tuple(parents.itervalues()) + (i,)
    
    @staticmethod
    def __sorted_trios(raw_trios):
        '''Find all genotyped trios in the genotype data set g. Use the pedigree information in p. The trios are
        sorted by parents.'''
        trios = np.array([trio for trio in raw_trios])
        return trios[np.lexsort((trios[:, 1], trios[:, 0]))] if trios.size else trios

    @staticmethod
    def __trios_to_families(trios):
        '''Collapse a collection of trios into a dictionary of parent-tuple-to-children-sets.'''
        family = {}
        for trio in trios:
            parents = tuple(trio[0:2])  # parents serve as family key
            family.setdefault(parents, im.Family(parents)).add(trio[2])  # add child to family children list
        return family
