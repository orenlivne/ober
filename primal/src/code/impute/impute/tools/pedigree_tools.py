'''
============================================================
Algorithms in pedigrees: queries, layout, kinship.

Created on May 31, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, tempfile, util, Image, shutil, famplot as fp, networkx as nx
from networkx.algorithms.traversal.depth_first_search import dfs_preorder_nodes
from numpy.ma.core import argmin
from numpy.core.numeric import Infinity
from impute.data.Person import Person
from impute.data.constants import MISSING, MATERNAL, PATERNAL
from collections import OrderedDict

#---------------------------------------------
# Constants
#---------------------------------------------
'''Maximum allowed generation gap between a parent and a child. This was empirically
found to be less than the minimum attainable value for the Hutterites pedigree, but
still serves as a correct guide from a human life expectancy viewpoint.'''
_MAX_GENERATION_GAP = 3

#---------------------------------------------
# Methods
#---------------------------------------------
def max_generation_gap(g, depth):
    '''Maximum generation gap between a parent and a child in the generation order depth.'''
    return max(depth[child] - depth[parent] for (parent, child) in g.edges_iter())

def node_depth(g, max_generation_gap=_MAX_GENERATION_GAP):
    '''Generate a dictionary of node depths for the pedigree p. 0-based. This is a topological
    order in the sense that parent depth < child depth for all nodes. 
    
    We also try to satisfy parent depth >= child depth - max_generation_gap by appropriately
    incrementing the depth of alien nodes (people that married into the pedigree but have no
    parent information so look like founder nodes); this feature is disabled if
    max_generation_gap is set to 0.'''
    
    # Create working copy of the graph
    h = g.copy()
    
    # Repeatedly set depth on 0-in-degree nodes, remove them, and increment depth
    d = 0
    depth = {}
    while (h.number_of_nodes() > 0):
        current_founders = [node for (node, in_degree) in h.in_degree_iter() if in_degree == 0]
        for node in current_founders:
            depth[node] = d
        h.remove_nodes_from(current_founders)
        d += 1

    # Fix alien nodes whenever possible        
    if max_generation_gap > 0: depth = __increment_alien_nodes(g, depth, max_generation_gap)
    
    return depth

def print_depth_report(g, depth=None):
    '''Print the list of nodes in the pedigree p, their depth and their parents' depths using the
    depth. If depth is not specified, the default is PedigreeTools.depth(p)'''
    if (depth == None):
        depth = node_depth(g)
    for child in g.nodes_iter():
        print 'Node %d: depth %d' % (child, depth[child]) 
        for (x, y) in ((parent, depth[parent]) for parent in g.predecessors_iter(child)):
            print '\tParent Node %d: depth %d' % (x, y)

def lowest_common_ancestor(g, u, v):
    '''Return (w,d), where a lowest common ancestor w of u and v with a minimum 
    d := dist(u,w)+dist(v,w), where dist(u,v) is the length of the path from u to v,
    and that distance d. If such a node is not found, returns (None, Infinity).'''
    # Initialize search at u and v
    d = 0
    u_parents = [u]
    v_parents = [v]
    # Keep track of distance of nodes from u and v
    U = dict(zip(u_parents, [d]))
    V = dict(zip(v_parents, [d]))

    # While not found an LCA, expand search to predecessor nodes 
    u_ancestors = set(u_parents)
    v_ancestors = set(v_parents)
    lca = u_ancestors & v_ancestors
    while (not lca and (u_parents or v_parents)):
        d += 1
        u_parents = all_predecessors(g, u_parents)
        v_parents = all_predecessors(g, v_parents)
        u_ancestors |= u_parents
        v_ancestors |= v_parents
        for node in u_parents: U[node] = d
        for node in v_parents: V[node] = d
        lca = u_ancestors & v_ancestors
    lca = list(lca)
    
    # Unique LCA or not found ==> finish
    if not lca: return None, Infinity
    elif (len(lca) == 1):
        w = lca[0]
        return w, U[w] + V[w]
    
    # Tie breaker of equal minimum distance sums.
    # If multiple such nodes are found, arbitrarily choose one of them.
    w = lca[argmin([U[node] + V[node] for node in lca])]
    return (w, U[w] + V[w]) if w else (None, Infinity)

'''Return the set of predecessor nodes of the node set nodeset in the DAG g.'''
all_predecessors = lambda g, nodeset: reduce(set.union, (set(g.predecessors(node)) for node in nodeset), set([]))

'''Return the set of successor nodes of the node set nodeset in the DAG g.'''
all_successors = lambda g, nodeset: reduce(set.union, (set(g.successors(node)) for node in nodeset), set([]))

def surrogate_parents(g, node, max_depth, min_depth=1, successors=True):
    '''Return a dictionary of surrogate-parent-to-depth of the node 'node' in the DAG g. Parents with
    min_depth <= depth <= max_depth are returned. If successors=False, successors of 'node' are omitted. '''
    # Build predecessor list. No need for a set here, since g is a DAG and predecessors of depth l
    # cannot intersect predecessors of depth < l
    result = OrderedDict()
    ancestors = [node]
    result[node] = 0
    for l in xrange(1, max_depth + 1):
        ancestors = reduce(list.__add__, (g.predecessors(x) for x in ancestors))
        result.update((x, l) for x in ancestors)
        if not ancestors: break
    
    # For each predecessor of depth l, find all successors of depth <= d-l. These can intersect.
    # m is the shortest path length between node and y, so do not update y's entry if it's already
    # in result with a lower m-value.
    for x, l in result.items():
        result.update((y, m) for (y, m) in ((y, l + len(path) - 1) for y, path in 
                                            nx.single_source_shortest_path(g, x, cutoff=max_depth - l).iteritems()
                                            if (True if successors else (node not in path)))
                                            if not result.has_key(y) or result[y] > m)
    # Remove entries of nodes with depth <= min_depth
    for y in (y for y, m in result.iteritems() if m < min_depth): del(result[y])
    return result

def families(g):
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

def find_families_by_either_parent(problem, parent, genotyped=True, min_children=0):
    '''find_families_by_parent(p, parent, min_children=0)
    
    Backward compability with older Problem NPZ archives that don\'t contain this method yet.
    Find families by parent ID.'''
    return [f for f in problem.families(genotyped=genotyped, min_children=min_children) 
            if parent in f.parents]
        
def family_member_union(families): 
    '''Return the union of all members in a family collection.'''
    try:
        return np.array(list(reduce(set.union, (f.member_set for f in families)))) 
    except TypeError:
        # Input to reduce() was empty
        return np.empty((0,), dtype=int)
    
def sibs(pedigree, node, parent_type):
    '''Return the indices of siblings (or half-siblings) of node 'node' on his parent_type-parent's side
    including i itself.'''
    parent = pedigree.parent(node, parent_type=parent_type, nullsafe=True)
    return pedigree.graph.successors(parent) if parent else []

def to_pedfiddler(pedigree, out_file, add_dummy=False, identifier='id', labels=None):
    '''Write pedigree to out file in pedfiddler input list format. If add_dummy=True,
    add dummy nodes for single-parent child missing parents. Sample identifier can be ''id'' or
    ''index''.'''
    if isinstance(out_file, str):
        with open(out_file, 'wb') as out:
            __to_pedfiddler(pedigree, out, identifier=identifier, labels=labels)            
    else:
        __to_pedfiddler(pedigree, out_file, identifier=identifier, labels=labels)

def draw_pedigree(pedigree, out_file, plot_params=fp.PlotParams(), colors=fp.DEFAULT_COLORS,
                  add_dummy=False, identifier='id', labels=None):
    '''Draw pedigree to file (EPS/PNG/...).'''
    extension = util.file_ext(out_file)
    if not extension in ['eps', 'png']:
        raise ValueError('Unsupported file extension ''%s''' % (extension,))

    coord_params = fp.CoordParams()
    coord_params.algorithm = 'default'
    
    # Create famplot input file
    # TODO: replace writing to a temp file by a copy constructor fp.Pedigree(our Pedigree) and
    # pass it instead of prefix
    tempdir = tempfile.mkdtemp()
    famplot_file = tempdir + '/file.fpl'
    eps_file = out_file if extension == 'eps' else tempdir + '/file.eps'

    to_pedfiddler(pedigree, famplot_file, add_dummy=add_dummy, identifier=identifier, labels=labels)
    p = fp.draw_pedigree(famplot_file, coord_params, plot_params, colors, eps_file)
    
    # Convert EPS to PNG/other supported formats by the Image library
    if extension != 'eps':
        with open(eps_file, 'rb') as eps:
            img = Image.open(eps)
            img.save(out_file)

    # Clean up
    shutil.rmtree(tempdir)
    return p

def draw_member_sub_pedigree(problem, i, out_file,
                             genotyped=False, plot_params=fp.PlotParams(), colors=fp.DEFAULT_COLORS,
                             add_dummy=False, identifier='id', labels=None):
    '''Draw the sub-pedigree comprising of all families that member i belongs to.'''      
    pedigree = problem.sub_pedigree(family_member_union(problem.find_families_by_member(i, genotyped=genotyped)))
    draw_pedigree(pedigree, out_file,
                  plot_params=plot_params, colors=colors, add_dummy=add_dummy, identifier=identifier, labels=labels)
    
def draw_member_neighbor_genotyped_pedigree(problem, i, depth, out_file,
                             plot_params=fp.PlotParams(), colors=fp.DEFAULT_COLORS,
                             add_dummy=False, identifier='id', labels=None):
    '''Draw the sub-pedigree comprising of i''s genotyped neighbors.'''
    nbhrs = problem.pedigree.neighbors_genotyped_selective(i, depth)
    # Add spouses of all quasi-founders to prevent multiple dummy nodes per spouse node in the pedigree drawing
    nbhrs = np.array(list(set([x for f in set([f for i in (x for x in nbhrs if not set(problem.pedigree.graph.predecessors(x)) & set(nbhrs)) for f in problem.pedigree.find_families_by_father(i, genotyped=False) + problem.pedigree.find_families_by_mother(i, genotyped=False)]) for x in f.member_list])))
    pedigree = problem.pedigree.sub_pedigree(nbhrs)
    draw_pedigree(pedigree, out_file,
                  plot_params=plot_params, colors=colors, add_dummy=add_dummy, identifier=identifier, labels=labels)
    return pedigree

'''Return the trios of the kids in the set, params.selected_samples, if in selected mode, otherwise,
all trios.'''
selected_trios = lambda problem, params: problem.kids_trios(params.selected_samples) if params.selected_mode else problem.trios()

'''Return the child-parent duos of the kids in the set, params.selected_samples, if in selected mode, otherwise,
all duos.'''
selected_duos = lambda problem, params, parent_type: problem.kids_duos(params.selected_samples, parent_type) if params.selected_mode else problem.duos(parent_type)

'''Return all families containing members of param.selected_samples, if in selected mode, else all
families.'''
selected_families = lambda problem, params, genotyped=False, min_children=0: \
    problem.find_families_by_members(params.selected_samples, genotyped=genotyped, min_children=min_children) \
    if params.selected_mode else problem.families(genotyped=genotyped, min_children=min_children)
    
#---------------------------------------------
# Private Methods
#---------------------------------------------
def __to_pedfiddler(pedigree, out, identifier='id', labels=None):
    '''Write pedigree to the file stream out in pedfiddler input list format. If add_dummy=True,
    add dummy nodes for single-parent child missing parents. Sample identifier can be ''id'' or
    ''index'' or ''node''.'''
    dummy_id = max(pedigree.sample_id)
    g = pedigree.graph
    
    # Add dummy parent nodes
    persons = [pedigree.person_with_ids(node) for node in g.nodes_iter()]
    persons_temp = dict(((pedigree.sample_id[node], pedigree.person_with_ids(node)) for node in g.nodes_iter()))
    for person in persons:
        f, m = person.father, person.mother
#        if (f == MISSING) and (m == MISSING):
#            # Quasi founder, ignore if it has no descendants, i.e., it is a disconnected node
#            # print person, pedigree.node_of[person.sample_id], pedigree.graph.successors(pedigree.node_of[person.sample_id])
#            if not pedigree.graph.successors(pedigree.node_of[person.sample_id]):
#                continue
            
        if (f == MISSING) ^ (m == MISSING):
            # If node has a single parent, create dummy node for the other parent 
            dummy_id += 1
            dummy_parent = Person(sample_id=dummy_id,
                                  sex=Person.SEX.MALE if (f == MISSING) else Person.SEX.FEMALE)
            persons_temp[dummy_id] = dummy_parent
            sibs = all_successors(g, g.predecessors(node))            
            # print 'dummy_parent', dummy_id, 'child', person.sample_id
            if f == MISSING:
                persons_temp[person.sample_id].father = dummy_id
                # print 'Updating child', person
                for child in sibs:
                    # print 'Updating sib', child, ' ', persons[child].sample_id
                    persons_temp[pedigree.sample_id[child]].father = dummy_id
            else:
                persons_temp[person.sample_id].mother = dummy_id
                # print 'Updating child', person
                for child in sibs:
                    persons_temp[pedigree.sample_id[child]].mother = dummy_id
        
    # Write pedigree to file
    for person in persons_temp.itervalues():
        __write_pedfiddler(pedigree, persons_temp, person, out, identifier, labels)

def __write_pedfiddler(pedigree, persons, person, out, identifier, labels):
    '''Write node line t out file in pedfiddler input list format.'''
    person_id = __person_identifier(pedigree, person, identifier)
    father_id = MISSING if person.father == MISSING else __person_identifier(pedigree, persons[person.father], identifier)
    mother_id = MISSING if person.mother == MISSING else __person_identifier(pedigree, persons[person.mother], identifier)
    # Custom identifier (dictionary of sample_id -> label
    label = (labels[person.sample_id] if labels.has_key(person.sample_id) else '-') if labels else str(person.sample_id)
    out.write("%7d %7d %7d %d %d %s\n" % (person_id, father_id, mother_id,
                                          person.sex, person.node_type, label))

def __person_identifier(pedigree, person, identifier):
    '''Convert a person''s sample_id to person''s identifier.'''
    if person.node_type == Person.TYPE.DUMMY:
        return person.sample_id
    elif identifier == 'node':
        return pedigree.node_of[person.sample_id]
    elif identifier == 'index':
        return pedigree.sample_index[pedigree.node_of[person.sample_id]]
    elif identifier == 'id':
        return person.sample_id
    else: 
        raise ValueError('Unsupported person identifier type ''%s''' % (identifier,))

def __increment_alien_nodes(g, depth, max_generation_gap):
    '''Try to satisfy parent depth >= child depth - max_generation_gap by appropriately
    incrementing the depth of alien nodes (people that married into the pedigree but have no
    parent information so look like founder nodes); this feature is disabled if
    max_generation_gap is set to 0.'''            

    # Identify pairs: alien node (people with no parent info that married into the
    # pedigree) + the other parent of a child they had together
    def aliens():
        for node in g.nodes_iter():
            if (g.out_degree(node) > 0):
                children = g.successors(node)
                # Node is more than max_generation_gap apart from all of its children
                if (depth[node] < min(depth[child] for child in children) - max_generation_gap):
                    # If at least one of the children has another parent from which we can infer 
                    # the proper depth, return the alien and that parent
                    for child in children:
                        if g.in_degree(child) > 1:
                            for parent in g.predecessors(child):
                                if (parent != node):
                                    yield (node, parent)
             
    # Decrement the depths of each alien and its ancestors to match the other parent's
    # depth whenever possible
    h = g.reverse(copy=True)
    
    for (alien, parent) in aliens():
        # We have information to update the alien node; see if we can actually increment its depth
        increment = depth[parent] - depth[alien]
        if (increment > 0):
            # Visit ancestors in depth-first-search pre-ordering (each ancestors is
            # processed before its parent)
            for ancestor in dfs_preorder_nodes(h, alien):
                new_depth = depth[ancestor] + increment
                if (g.out_degree(ancestor) == 0) or (new_depth < min(depth[node] for node in g.successors(ancestor))):
                    # Not violating topological order ==> update this ancestor
                    depth[ancestor] = new_depth
    
    return depth

__PARENT_TYPE_STR = {PATERNAL: 'PATERNAL', MATERNAL: 'MATERNAL'}

def parent_type_str(parent_type):
    '''Pretty-print a parent type constant.'''
    return __PARENT_TYPE_STR[parent_type] if parent_type else 'NONE'

