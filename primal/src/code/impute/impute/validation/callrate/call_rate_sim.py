#!/usr/bin/env python
'''
============================================================
Perform gene dropping simulation to estimate the ideal full-
and partial-genotype call rates from the pedigree.

Prerequisite dependencies:
* networkx - http://networkx.github.com/download.html
* numpy - http://www.numpy.org/

Input: pedigree text file
Row format: <id> <father_id> <mother_id> <type>
type=0 for samples to be ignored in the output
type=1 for samples to estimate call rates for
type=2 for typed samples (e.g., the 98 WGS Hutterites samples) 
id=0 indicates missing information (e.g., founders' parents)

Created on March 12, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
from __future__ import division
import itertools, os, sys, optparse, csv, networkx as nx, numpy as np, random
from networkx.algorithms.dag import is_directed_acyclic_graph

#---------------------------------------------
# Constants
#---------------------------------------------
'''Indicates missing allele data in a person.'''
INVALID_ALLELE = -1

'''Indicates missing data in adjacency lists.'''
MISSING = 0

'''Paternally-inherited allele index.'''
PATERNAL = 0
'''Maternally-inherited allele index.'''
MATERNAL = 1
'''Both allele locations'''
ALLELES = (PATERNAL, MATERNAL)

'''Sample types. Target=to be imputed. Typed=WGS samples.'''
class NODE_TYPE: IGNORED, TARGET, TYPED = xrange(3)

#---------------------------------------------
# Methods
#---------------------------------------------
def call_rate_sim(options):
    '''Main call rate simulation program.'''
    # Load pedigree and ID data
    ped = Pedigree.load(options.ped_file)
    
    # Run simulations
    sims = [ImputationSimulation(ped, debug=options.debug, alpha=alpha) for alpha in linspace(0.9, 1.0, options.num_alphas)]
    for sim in sims:
        if options.debug:
            sys.stdout.write('Simulation, phasing %% %.3f\n' % (sim.alpha))
        sim.run(options.num_simulations)
    return sims

def parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    prog = os.path.basename(sys.argv[0])
    usage = 'Usage: %s [flags] <ped_file>\n\n' \
    'Perform gene dropping simulation to estimate the ideal full-\n' \
    'and partial-genotype call rates from the pedigree.\n' \
    '\n' \
    'Prerequisite dependencies:\n' \
    '* networkx - http://networkx.github.com/download.html\n' \
    '\n' \
    'Input: pedigree text file\n' \
    'Row format: <id> <father_id> <mother_id> <type>\n' \
    'type=0 for samples to be ignored in the output\n' \
    'type=1 for samples to estimate call rates for\n' \
    'type=2 for typed samples (e.g., the 98 WGS Hutterites samples) \n' \
    'id=0 indicates missing information (e.g., founders'' parents)\n\n' \
    'Type ''%s -h'' to display full help.' % (prog, prog)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-n', '--num-simulations', type='int'  , dest='num_simulations',
                      default=100,
                      help='Number of simulations to run')
    parser.add_option('-a', '--num-alphas', type='int'  , dest='num_alphas',
                      default=5,
                      help='Number of phasing %% values between 0.9 and 1.0 to run for')
#    parser.add_option('-e', '--target-relative-error', type='float'  , dest='e',
#                      default=0.01,
#                      help='Target relative error in call rates')
    (options, args) = parser.parse_args(sys.argv[1:])
    
    # Argument validation
    if len(args) != 1:
        print usage
        sys.exit(1)
        
    # Read file name
    options.ped_file = args[0]
    return options

haps = lambda x: set(itertools.product(x, ALLELES))

def random_allele():
    '''Return the paternal or maternal allele index with equal probabilities.'''
#    return random.randint(0, 1)
    return np.random.randint(0, 2)

def linspace(start, stop, n):
    if n == 1:
        yield stop
        return
    h = (stop - start) / (n - 1)
    for i in xrange(n): yield start + h * i
        
####################################################################################
class Pedigree(object):
    '''
    ============================================================
    Pedigree data structure. This is a networkx directed graph
    with Person-type nodes.
    ============================================================
    '''
    def __init__(self, graph):
        '''Construct a pedigree from a Person DAG (edge goes from parent to child).'''
        
        # Validate graph
        if (not is_directed_acyclic_graph(graph)):
            raise ValueError('Pedigree graph must be a directed acyclic graph (DAG)')
        for node in graph.nodes_iter():
            degree = graph.in_degree(node)
            if (degree > 2):
                raise ValueError('In-degree must be at most 2, but degree[%d] = %d' % (node, degree,))
        
        # Transform data to desired form: add dummy founder nodes
        self.graph = graph
        self.graph = self.__add_dummy_founders()
        
        # Lazily-initialized cached properties
        self._sorted_nodes = None 
        self._founder_nodes = None
        self._child_nodes = None
        
        g = self.graph
        for person in self.child_nodes:
            f, m = self.father(person), self.mother(person)
            g.node[person]['father'] = g.node[f]['alleles']
            g.node[person]['mother'] = g.node[m]['alleles']
#            print person, self.father(person), self.mother(person)
#            print g.node[f]['alleles']
#            print g.node[person]['father']
#            g.node[f]['alleles'][0] = 1000
#            print g.node[f]['alleles']
#            print g.node[person]['father']
            
        # Generate founder alleles (0-based). These are fixed for all simulations.
        k = 0
        for person in self.founder_nodes:
            alleles = g.node[person]['alleles']
            alleles[PATERNAL] = k
            alleles[MATERNAL] = k + 1
            k += 2
        self.num_founder_alleles = k
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    @staticmethod
    def load(infile, separator=' '):
        '''Read pedigree from text file.
        Expected row format: <id> <father_id> <mother_id> <type>'''
        # Build nodes=samples and edges=parent-to-offspring; indicate which parent on the edge
        g = nx.DiGraph()
        with open(infile, 'rb') as infile:
            for line in csv.reader(infile, delimiter=separator, skipinitialspace=True):
                person, father, mother, node_type = map(int, line[0:4])
                # print person, father, mother, node_type
                g.add_node(person, dict(alleles=[INVALID_ALLELE, INVALID_ALLELE], node_type=node_type, father=None, mother=None))
                if father != MISSING:
                    g.add_edge(father, person, dict(parent=PATERNAL))
                if mother != MISSING:
                    g.add_edge(mother, person, dict(parent=MATERNAL))
        # print sorted(g.nodes())
        return Pedigree(g)
    
    def __repr__(self):
        '''Succinct textual representation of the pedigree.'''
        return 'Pedigree[nodes=%d, edges=%d]' % (self.graph.number_of_nodes(), self.graph.number_of_edges()) 

    def pprint(self):
        '''Pretty-print the pedigree.'''
        print 'Pedigree[nodes=%d, edges=%d]' % (self.graph.number_of_nodes(), self.graph.number_of_edges())
        print 'Node F  M  Type'
        for (name, person) in self.graph.nodes_iter():
            print '%6d %6d %6d' % (name, self.father(person), self.mother(person), person.type) 

    def nodes_of_type(self, node_type):
        '''Return all nodes of type node_type.'''
        return list(filter(lambda x: self.graph.node[x]['node_type'] == node_type, self.graph.nodes_iter()))
    
    def father(self, person):
        '''Father node of node ''person''.'''
        return self.parent(person, PATERNAL)

    def mother(self, person):
        '''Mother node of node ''person''.'''
        return self.parent(person, MATERNAL)

    def parent(self, person, parent_type):
        '''Parent node.'''
        parent = [e[0] for e in self.graph.in_edges(person, True) if e[2]['parent'] == parent_type]
        return parent[0] if parent else MISSING

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def sorted_nodes(self):
        '''Return the nodes in topological order (parent precedes child).'''
        if self._sorted_nodes is None:
            self._sorted_nodes = nx.topological_sort(self.graph)
        return self._sorted_nodes

    @property
    def founder_nodes(self):
        '''Return the founder nodes (nodes without parents).'''
        if self._founder_nodes is None:
            self._founder_nodes = filter(lambda x: self.graph.in_degree(x) == 0, self.graph.nodes_iter())
        return self._founder_nodes

    @property
    def child_nodes(self):
        '''Return the non-founder nodes (nodes with parents), in topological order.'''
        if self._child_nodes is None:
            self._child_nodes = filter(lambda x: self.graph.in_degree(x) > 0, self.sorted_nodes)
        return self._child_nodes
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __add_dummy_founders(self):
        '''Add dummy founders to children with a single parent.'''
        # Need a copy, since we would be changing g during iterating over g's nodes, which is a no-no
        g = self.graph
        g_new = g.copy()
        # Increment dummy node names starting at the (largest ID in the current pedigree)+1
        dummy_founder = max(g.nodes_iter())
        
        for person in g.nodes_iter():
            # If node has a single parent, create dummy node for the other parent
            # print person, self.father(person), self.mother(person)
            # print person, self.graph.in_degree(person)
            if self.graph.in_degree(person) == 1:
                dummy_founder += 1
                f = self.father(person)
                # print person, f, m 
                parent_type = PATERNAL if (f == MISSING) else MATERNAL
                # Add dummy parent - one per each node, since we don't have information that
                # they share the same missing parent 
                g_new.add_node(dummy_founder, dict(alleles=[INVALID_ALLELE, INVALID_ALLELE], node_type=NODE_TYPE.IGNORED))
                g_new.add_edge(dummy_founder, person, dict(parent=parent_type))
        return g_new

####################################################################################
class ImputationSimulation(object):
    '''Simulates gene dropping and measures call rates (based on a pedigree)
    between a target sample ID i and each of the WGS training samples j in T.'''

    #---------------------------------------------
    # Constants
    #---------------------------------------------
    # Paternal allele = left bit; maternal allele = right bit in the 'imputable' dictionary values below
    _ALLELE_MASK = {PATERNAL: 2, MATERNAL: 1}

    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, ped, alpha=1.0, debug=False):
        '''e=Desired relative error in call rate estimates.'''
        # Pedigree
        self.ped = ped
        self.typed_haps = haps(ped.nodes_of_type(NODE_TYPE.TYPED))
        self.target = ped.nodes_of_type(NODE_TYPE.TARGET)
        self.target_haps = haps(self.target)
        self.alpha = alpha
        self.debug = debug
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def init(self, num_simulations):
        '''Reset the simulation.'''
        # Genotype call counter for each target sample (a tuple of (count of fully-called genotypes, count of genotypes with one allele called)
        self.count_sample = dict((person, [0, 0]) for person in self.target) 
        # Call rates per SNP
        self.count_snp = [[0, 0] for _ in xrange(num_simulations)]
        # Total # simulations
        self.count = 0

    def run(self, num_simulations):
        '''Run n simulations. If n=None, uses default #.'''
        self.init(num_simulations)
        m = max(num_simulations / 10, 1000)
        for i in xrange(num_simulations):
            self._simulate_imputation()
            if self.debug and i % m == 0:
                print '\t%d/%d simulations completed ...' % (i, num_simulations)

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def call_rate_sample(self):
        '''Full genotype call rate, partial gentype call rate (%), allele call rate estimates
        for each target sample.'''
        return dict((k, ((1.0 * v[0]) / self.count,
                         (1.0 * (v[0] + v[1])) / self.count,
                         (2.0 * v[0] + v[1]) / (2 * self.count)))
                     for k, v in self.count_sample.iteritems())
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------        
    def _simulate_imputation(self):
        '''Simulate imputing i from phased T-haplotypes.'''
        target, target_haps = self._target_and_target_haps()
        
        # Simulate gene dropping from parent to child
        ped, g = self.ped, self.ped.graph
        for person in ped.child_nodes:
            # print person, g.node[person]
            alleles = g.node[person]['alleles']
            alleles[PATERNAL] = g.node[person]['father'][random_allele()]
            alleles[MATERNAL] = g.node[person]['mother'][random_allele()]

        # Assign haplotypes to IBD sets
        ibd = [set([]) for _ in xrange(self.ped.num_founder_alleles)]
        for person in ped.sorted_nodes:
            # print person, g.node[person]['alleles']
            for parent_type, allele in enumerate(g.node[person]['alleles']):
                # print 'ibd[%d] <- (%d,%d)' % (allele, person, parent_type)
                ibd[allele].add((person, parent_type))
        
        # Infer which alleles of the target samples can be imputed from IBD-sharing typed samples
        imputable = dict((k, 0) for k in target)
        for person, allele in (x for ibd_set in ibd if ibd_set & self.typed_haps for x in ibd_set if x in target_haps):
            imputable[person] |= ImputationSimulation._ALLELE_MASK[allele]
        
        # Update counts
        self.count += 1
        # Both alleles are imputable
        for i, person in enumerate(person for person, value in imputable.iteritems() if value == 3):
            self.count_sample[person][0] += 1
        self.count_snp[0] = i + 1
        
        # Either allele is imputable, but not both 
        for i, person in enumerate(person for person, value in imputable.iteritems() if value == 1 or value == 2):
            self.count_sample[person][1] += 1
        self.count_snp[1] = i + 1

    def _target_and_target_haps(self):
        '''Return the target sample list and target haplotype set for a simulation. Delete (1-alpha) of
        the individuals in the original target set.'''
        if abs(self.alpha - 1.0) < 1e-15:
            return self.target, self.target_haps
        else:
            restricted_target = random.sample(self.target, int(round(self.alpha * len(self.target))))
            return restricted_target, haps(restricted_target)

####################################################################################
if __name__ == '__main__':

    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    options = parse_command_line_args()
    sims = call_rate_sim(options)

    # Print results
    for person in sims[0].call_rate_sample.iterkeys():
        sys.stdout.write('%-10d' % (person,))
        for sim in sims:
            call_rate = sim.call_rate_sample[person]
            sys.stdout.write(' %f %f' % (call_rate[0], call_rate[2]))
        sys.stdout.write('\n')
