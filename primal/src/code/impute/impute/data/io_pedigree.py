'''
============================================================
Read and write Pedigree objects to/from PLINK TFAM files.

Created on May 31, 2012
Rerganized on August 4, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, networkx as nx, util
from impute.data.Pedigree import Pedigree
from impute.data.constants import MISSING, PATERNAL, MATERNAL
from impute.data.Person import Person
from impute.data import constants

#---------------------------------------------
# Constants
#---------------------------------------------
'''See plink format documentation at http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
Dummy hard-coded values for now.'''
_FAMILY_ID = 'HUTTERITES'

#---------------------------------------------
# Methods
#---------------------------------------------
def read(file_name, genotyped_id_file=None):
    '''Load a pedigree from file.'''
    # Read optional genotyped sample ids
    genotyped_ids = read_sample_id(genotyped_id_file) if genotyped_id_file is not None else None
    return __read_pedigree(file_name, genotyped_ids=genotyped_ids)

'''Read sample IDs from a TFAM file into a numpy array.'''
read_sample_id = lambda tfam: np.genfromtxt(tfam, dtype=np.int)[:, 1]

'''Return an array whose ith element is the genotyped_pedigree sample index of the ith sample index in sub_pedigree.'''
genotyped_sample_index = lambda genotyped_pedigree, sub_pedigree: np.array(map(genotyped_pedigree.node_of.get, sub_pedigree.sample_id[0:sub_pedigree.num_genotyped]))

####################################################################################
def write(pedigree, file_name, save_node_type=True):
    '''Save pedigree to file. file_name is a file set prefix or a file descriptor
    Note: sample_ids are not saved if file_name is a file descriptor.'''
    if isinstance(file_name, str):
        # file_name is a file set prefix
        with open(file_name + '.pdg.tfam', 'wb') as out_file:
            __write_pedigree(pedigree, out_file)
        with open(file_name + '.tfam', 'wb') as out_file:
            __write_pedigree(pedigree, out_file, genotyped_only=True, save_node_type=save_node_type)
    else:
        # file name is a file descriptor
        __write_pedigree(pedigree, file_name, save_node_type=save_node_type)
        
#---------------------------------------------
# Private Methods
#---------------------------------------------
def __read_pedigree(file_name, genotyped_ids=None):
    '''Load pedigree from the reads from the PLINK TFAM file file_name. If a list of genotyped ids
    genotyped_ids is specified, the pedigree node IDs are recoded to consecutive non-negative
    integers, i.e., the original plink ID list is mapped to the node list 1..(len(plink_id)).
    
    If genotype IDs are available, set node type to genotyped/not genotyped using those values;
    otherwise, fall back to the data[:,5] column.''' 
    
    # Load data from text file_name
    data = np.genfromtxt(file_name, dtype=np.dtype(int), usecols=range(1, 6))
    nodes, missing = data[:, 0], MISSING
    num_samples, num_columns = data.shape
    if genotyped_ids is None:
        sample_id = nodes
        order = np.arange(0, num_samples)
        num_genotyped = len(sample_id)
        node_type = data[:, 5] if num_columns >= 6 else np.tile(constants.INDETERMINATE, (num_samples, 1))
        node_to_sample_id = None
    else:
        # Recode nodes if a genotyped ID list was specified

        # Dummy node ID for missing data in the pedigree adjacency list. Must not be a
        # possible node number (typically, a negative number should work).'''
        missing = -1
        if isinstance(genotyped_ids, np.ndarray):
            genotyped_ids = genotyped_ids.tolist()
        rest = sorted(set(nodes) - set(genotyped_ids))
        # Keep track of original IDs
        sample_id = np.array([missing] + genotyped_ids + list(rest))
        recoding = dict(zip(sample_id, [missing] + range(Pedigree.START_ID, Pedigree.START_ID + data.shape[0])))
        # Recode ID columns (first three); since array is thin and long, recode columns and
        # then transpose
        num_id_columns = 3
        data = np.concatenate((np.array([[recoding[x] if recoding.has_key(x) else missing for x in data[:, j]]
                                         for j in xrange(0, num_id_columns)]).transpose(),
                               data[:, num_id_columns:]), axis=1)
        nodes = data[:, 0]
        # Remove the sample_id entry of the missing value 
        # sample_id = sample_id[1:] 
        node_to_sample_id = util.dict_invert(recoding)
        
        num_genotyped = len(genotyped_ids)
        node_type = np.array([(Person.TYPE.GENOTYPED if x < num_genotyped else Person.TYPE.NOT_GENOTYPED)
                     for x in xrange(0, num_samples)])

    # Construct _graph
    graph = nx.DiGraph()
    # Add all nodes
    graph.add_nodes_from(nodes)
    # Add father->child edges, mother->child edges
    for (parent_type, column) in {PATERNAL: 1, MATERNAL: 2}.iteritems():
        graph.add_edges_from(data[:, (column, 0)], type=parent_type)
    # Remove missing data = edges from nodes to the dummy node 'missing'
    graph.remove_node(missing)
    
    if genotyped_ids is not None:
        # After graph might have reordered nodes, reorder sample_ids accordingly
        node_to_line = dict(zip(nodes, graph.nodes()))
        sample_id = np.array([node_to_sample_id[x] for x in graph.nodes()])
        order = np.array([node_to_line[x] for x in graph.nodes()])

    # Check that in-degree is at most 2
    if np.nonzero(np.array(graph.in_degree().values()) > 2)[0]:
        raise ValueError('Input data contains a child node with more than two parents')

    # If genotype IDs are available, set node type to genotyped/not genotyped using those values;
    # otherwise, fall back to the data column
    return Pedigree(graph, sample_id=sample_id, num_genotyped=num_genotyped,
                    sex=data[order, 3], phenotype=data[order, 4], node_type=node_type)

####################################################################################
def __write_pedigree(pedigree, out_file, genotyped_only=False, save_node_type=True):
    '''Write a pedigree to PLINK TFAM file. If genotyped_only=True, output only genotyped_nodes
    and preserve their original order. Otherwise, output all nodes, sorted by ID.'''
    fmt = __line_format(save_node_type)
    nodes = pedigree.graph.nodes()
    if genotyped_only:
        lines = (__output_line(pedigree, i, save_node_type) for i in nodes[0:pedigree.num_genotyped])
    else:
        lines = sorted((__output_line(pedigree, i, save_node_type) for i in nodes),
                        key=lambda x: x[1])
    for line in lines:
        out_file.write((fmt + '\n') % line)
            
def __output_line(pedigree, i, save_node_type):
    '''Generates output line of node i.'''
    person = pedigree.person_with_ids(i)    
    fields = (_FAMILY_ID, person.sample_id, person.father, person.mother, person.sex, person.phenotype)
    if save_node_type:
        fields += (person.node_type,)
    return fields

def __line_format(save_node_type):
    '''Line output format.'''
    fmt = '%s %d %d %d %d %d'
    if save_node_type:
        fmt += ' %d'
    return fmt
