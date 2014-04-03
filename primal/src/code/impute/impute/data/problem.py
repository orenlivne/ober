'''
============================================================
A genetic problem to solve. This is the request type passed
betweened phasing stages. it contains pedigree, genotype and
haplotype objects, as well as other attributes useful in the
course of computation.

Created on July 6, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, util, impute as im

####################################################################################
class Problem(object):
    '''The overall problem to solve for all persons and all snps.'''
    #---------------------------------------------
    # Constants
    #---------------------------------------------

    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, pedigree, genotype, haplotype=None, error=None, info=None,
                 frames=None, lam=None):
        '''Allocate a problem with a pedigree, genotype and empty haplotype set
        (or the specified initial haplotype).'''
        self.pedigree = pedigree
        self.genotype = genotype
        self.haplotype = haplotype if haplotype is not None else \
        im.factory.GenotypeFactory.empty_from_genotype(genotype)
        # An array where a non-zero indicates a flagged genotype error and may encode the source
        # of the error or the stage at which it was flagged
        self.error = error if error is not None else np.zeros((self.num_snps, self.num_samples), dtype=np.int8)

        # SNP Meta data
        self.frames = frames

        # Sample kinships
        self.lam = lam 

        self._info = ProblemInfo(pedigree, genotype) if info is None else info
        
        # Lazily-initialized cached properties        
        self._snp_range = None
        self._lam_eval = None
                
    #---------------------------------------------
    # Operators
    #---------------------------------------------
    def __key(self):
        '''Uniquely-identifying key of this object.'''
        return (self.pedigree, self.genotype, self.error.tolist(), self.frames, self._info)

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
        return 'Problem[snps=%d, samples=%d, pedigree=%d, errors=%d]' \
            % (self.num_snps, self.num_samples, self.num_pedigree_samples, self.num_errors)
    
    def __add__(self, other):
        '''Merge two SegmentComposites (binary operation).'''
        '''Merge the SNPs of two problems with the same sample set. Return the merged Problem object.
        Error, QC arrays and frames are omitted.'''        
        g = np.concatenate((self.g, other.g), axis=0)
        g_snp = np.concatenate((self.genotype.snp, other.genotype.snp))
        genotype = im.factory.GenotypeFactory.new_instance('genotype', g, g_snp)

        h = np.concatenate((self.h, other.h), axis=0)
        h_snp = np.concatenate((self.haplotype.snp, other.haplotype.snp))
        haplotype = im.factory.GenotypeFactory.new_instance('haplotype', h, h_snp)

        snp = np.concatenate((self.info.snp, other.info.snp))
        info = ProblemInfo(self.pedigree, genotype, snp=snp)
        return Problem(self.pedigree, genotype, haplotype, info=info, frames=None, lam=self.lam)
    
    def __radd__(self, other):
        '''Reverse-add a problem to this object (binary operation).'''
        if not other:
            # Empty + self = self
            return self
        elif isinstance(other, self.__class__):
            # SegmentComposite (of the same sub-class type as ours) is the only supported class
            return self +other
        else:
            raise ValueError('Cannot add an object of type SegmentComposite to a %s' % (other.__class__,))
     
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def genotype_error(self, snps, samples, msg):
        '''Flag a genotype error. If scalar arguments are passed in, they are regarded as a single error.'''  
        
        # if np.isscalar(snp):
        #    (snp, sample) = ([snp], [sample])
        # for (i, s) in enumerate(snp):
        #    self._errors.add((s, sample[i]))
        
        # log(WARN, msg)
        if snps.size:
            # print '| %-20s | SNP %s sample %s %s' % ('Error', repr(snps), repr(samples), msg)
            self.genotype.clear(snps, samples)
            self.haplotype.clear(snps, samples)
            self.error[snps, samples] = im.constants.ERROR

    ####################################################################################
    # Sample Kinships
    ####################################################################################
    def lam_eval(self, f):
        '''Return the recombination rate lambda(f) where f=inbreeding (or kinship) coefficient.'''
        if not self._lam_eval: self._lam_eval = LambdaEvaluator(self.lam)
        return self._lam_eval(f)

    ####################################################################################
    # Statistics
    ####################################################################################
    def is_genotyped(self, x):
        '''Does the sample index x correspond to a genotyped sample or not.''' 
        return self.pedigree.is_genotyped(x)
    
    def fill_fraction(self, snps=None, snp_range=None, sample=None, allele=None):
        '''Return an array whose columns are (sample index, filled haplotype %), for all sample
        indices (by default) or the selected set of sample indices.'''
        sample = ([sample] if np.isscalar(sample) else sample) if sample is not None else xrange(0, self.num_samples)
        fill = np.array([(x, self.haplotype.fill_fraction(snps=snps, snp_range=snp_range, sample=x, allele=allele)) 
                        for x in sample if self.is_genotyped(x)])
        # fill = np.array([(self.pedigree.sample_id[x], self.haplotype.fill_fraction(sample=x)) 
        #                for x in sample])
        # Sort by ascending fill %
        # fill = fill[np.argsort(fill[:,1])]
        return fill

    def fill_fraction_of_sample(self, sample):
        '''Single-sample hap fill fraction. (A useful alias.)'''
        return self.haplotype.fill_fraction(sample=sample)
    
    def find_samples_with_fill_ge(self, fill_threshold, sample=None):
        '''Return the sample IDs (column 0) and fill fraction (column 1) of all samples
        (or a subset 'sample' of the samples, if sample is not None) with fill fraction => fill_threshold.'''
        f = self.fill_fraction(sample=sample)
        return f[np.where(f[:, 1] >= fill_threshold)[0], :] if f.size else np.empty((0, 1)) 

    def find_samples_with_fill_lt(self, fill_threshold, sample=None):
        '''Return the sample IDs (column 0) and fill fraction (column 1) of all samples
        (or a subset 'sample' of the samples, if sample is not None) with fill fraction => fill_threshold.'''
        f = self.fill_fraction(sample=sample)
        return f[np.where(f[:, 1] < fill_threshold)[0], :] if f.size else np.empty((0, 1)) 

    ####################################################################################
    # Duos, Trios & Family Lists
    ####################################################################################
    def duos(self, parent_type, genotyped=True):
        '''Find and cache all genotyped duos in the data set:
        (child-parent[parent_type]) tuples.'''
        return self.pedigree.duos(parent_type, genotyped=genotyped)
    
    def trios(self, genotyped=True):
        '''Find and cache all genotyped trios in the data set. Sorted by parents.'''
        return self.pedigree.trios(genotyped=genotyped)

    def kids_duos(self, kids, parent_type, genotyped=True):
        '''Find and cache all genotyped duo IDs in the data set:
        (child-parent[parent_type]) tuples where child is in kids. Not cached.'''
        return self.pedigree.kids_duos(kids, parent_type, genotyped=genotyped)

    def kids_trios(self, kids, genotyped=True):
        '''Return the trios of the kids in the list ''kids''. Not cached.'''
        return self.pedigree.kids_trios(kids, genotyped=genotyped)

    def families(self, genotyped=True, min_children=0, max_children=np.inf): 
        '''Return an iterator over families in a list of trios in which both parents are in G and at
        least min_children of their children are in G. min_children=0 by default. Cached.'''
        return self.pedigree.families(genotyped=genotyped, min_children=min_children, max_children=max_children)
    
    def families_union(self, genotyped=True, min_children=0, max_children=np.inf): 
        '''Return the union of all genotyped nuclear family members with at least min_children children
        in the trio list trio.'''
        return self.pedigree.families_union(genotyped=genotyped, min_children=min_children, max_children=max_children)
    
    ####################################################################################
    # Family Info and Selection
    ####################################################################################
    def set_family_info(self, family, info, genotyped=True):
        '''Family-to-info dictionary entry setting.'''
        self._info.set_family_info(family, info)

    def find_family(self, father, mother, genotyped=True):
        return self.pedigree.find_family(father, mother, genotyped=genotyped)
        
    def find_families_by_either_parent(self, parent, genotyped=True, min_children=0):
        '''Find families by parent ID.'''
        return self.pedigree.find_families_by_either_parent(parent, genotyped=genotyped, min_children=min_children)

    def find_families_by_parent(self, parent_type, parent, genotyped=True, min_children=0):
        '''Find families by parent type + ID.'''
        return self.pedigree.find_families_by_parent(parent_type, parent, genotyped=genotyped, min_children=min_children)
        
    def find_families_by_father(self, father_id, genotyped=True, min_children=0):
        '''Find families by father ID.'''
        return self.pedigree.find_families_by_father(father_id, genotyped=genotyped, min_children=min_children)
        
    def find_families_by_mother(self, mother_id, genotyped=True, min_children=0):
        '''Find families by mother ID.'''
        return self.pedigree.find_families_by_mother(mother_id, genotyped=genotyped, min_children=min_children)

    def find_family_by_child(self, child_id, genotyped=True, min_children=0):
        '''Find families by child ID.'''
        return self.pedigree.find_family_by_child(child_id, genotyped=genotyped, min_children=min_children)

    def find_families_by_member(self, member_id, genotyped=True, min_children=0):
        '''Find all families that contain a member ID.'''
        return self.pedigree.find_families_by_member(member_id, genotyped=genotyped, min_children=min_children)

    def find_families_by_members(self, member_ids, genotyped=True, min_children=0):
        '''Find all families that contain a member ID.'''
        return self.pedigree.find_families_by_members(member_ids, genotyped=genotyped, min_children=min_children)
    
    ####################################################################################
    # Sub-problems
    ####################################################################################
    def sub_problem_of_parents(self, father, mother, snps=None):
        f = self.find_family(father, mother, genotyped=False)
        if not f:
            raise ValueError('Family (%d,%d) not found' % (father, mother))
        return self.sub_problem_of_family(f, snps=snps)
     
    def sub_problem_of_family(self, f, snps=None):
        '''Return a sub-problem that contains a single nuclear genotyped family (father,mother)
        and a subset of the SNPs (or all SNPs, if snps=None). In the sub-problem,
        nodes are a running sequence and sample_ids refer to sample indices in THIS problem.'''
        problem = self.sub_problem(f.member_list, snps=snps)
        
        # Add in family information, if available
        if self.info.family_info.has_key(f):
            new_f = problem.families(genotyped=False)[0]  # Convert to from self to problem's coordinates
            self_info = self.info.family_info[f]
            family_info = {new_f: {im.constants.PATERNAL: None, im.constants.MATERNAL: None}}
            for parent_type in im.constants.ALLELES:
                entry = self_info[parent_type]
                if entry:
                    family_info[new_f][parent_type] = entry.copy()
            problem.info.family_info = family_info
        return problem

    def sub_pedigree(self, samples):
        '''Return the sub-pedigree containing the selected samples only.'''
        return self.pedigree.sub_pedigree(samples)
    
    def sub_problem_of_snps(self, snps=None):
        '''Return a sub-problem that contains a subset of the SNPs. Assumes at least two SNPs.'''
        # Re-order samples so that all genotyped appear before all non-genotyped

        # Create deep copies of relevant parts of data arrays
        g, h = self.data
        g_snp, h_snp, qc = self.genotype.snp, self.haplotype.snp, self.haplotype.qc
        frames = self.frames 
        if snps is not None:
            g, h = g[snps, :, :], h[snps, :, :]
            if qc.size: qc = qc[snps, :, :]
            g_snp, h_snp = g_snp[snps], g_snp[snps]
            # Restrict frames to snps, convert to new SNP indices
            orig_snp = dict((v, k) for k, v in enumerate(snps))
            def sub_frame(frame):
                for x in frame:
                    if orig_snp.has_key(x): yield orig_snp[x]
            
            frames = util.mdict()
            for k, v in frames.iteritems():
                for frame in v:                    
                    frames[k] = sub_frame(frame)
            
        g, h = g.copy(), h.copy()
        if qc.size:
            qc = qc.copy()
        g_snp, h_snp = g_snp.copy(), h_snp.copy()
        
        # Build sub-problem object graph
        g = im.factory.GenotypeFactory.new_instance('genotype', g, g_snp)
        h = im.factory.GenotypeFactory.new_instance('haplotype', h, h_snp)
        h.qc = qc
        if self.haplotype.poo_phase is not None:
            h.poo_phase = self.haplotype.poo_phase.copy()
        if self.haplotype.hap_type is not None:
            h.hap_type = self.haplotype.hap_type[snps].copy()
        
        # Build restricted info object
        error = self.error.copy() if self.error.size else self.error
        if error.size: error = error[snps, :]
        ibd = im.segment.SegmentSet(im.segment.Segment(x.snp, map(lambda y: (self.sample_index[y[0]], y[1]), x.samples),
                                                       x.bp, error_snps=x.error_snps) for x in self.info.ibd)
        info = ProblemInfo(self.pedigree, g, snp=(self.info.snp[snps] if snps is not None else self.info.snp), ibd=ibd)
        return Problem(self.pedigree, g, haplotype=h, info=info, error=error, frames=frames, lam=self.lam)

    def sub_problem(self, samples, snps=None):
        '''Return a sub-problem that contains a subset 'samples' of the genotyped nodes.'''
        # Re-order samples so that all genotyped appear before all non-genotyped
        if isinstance(samples, list): samples = np.array(samples)
        # Create pedigree object
        p = self.sub_pedigree(samples)

        genotyped = np.where(samples < self.pedigree.num_genotyped)[0]
        num_genotyped = len(genotyped) 
        samples = samples[np.concatenate((genotyped, np.where(samples >= self.pedigree.num_genotyped)[0]))]

        # Create deep copies of relevant parts of data arrays
        g, h = self.data
        g_snp, h_snp, qc = self.genotype.snp, self.haplotype.snp, self.haplotype.qc
        frames = self.frames 
        if snps is not None:
            g, h, qc = g[snps, :, :], h[snps, :, :], qc[snps, :, :] if qc.size else None
            g_snp, h_snp = g_snp[snps], g_snp[snps]
            # Restrict frames to snps, convert to new SNP indices
            orig_snp = dict((v, k) for k, v in enumerate(snps))
            def sub_frame(frame):
                for x in frame:
                    if orig_snp.has_key(x): yield orig_snp[x]
            
            frames = util.mdict()
            for k, v in frames.iteritems():
                for frame in v:                    
                    frames[k] = sub_frame(frame)
            
        genotyped = p.sample_index[0:num_genotyped]
        g, h = g[:, genotyped, :].copy(), h[:, genotyped, :].copy()
        if qc.size: qc = qc[:, genotyped, ].copy()
        g_snp, h_snp = g_snp.copy(), h_snp.copy()
        
        # Build sub-problem object graph
        g = im.factory.GenotypeFactory.new_instance('genotype', g, g_snp)
        h = im.factory.GenotypeFactory.new_instance('haplotype', h, h_snp)
        h.qc = qc
        # Build restricted info object
        error = self.error[:, genotyped].copy() if self.error.size else self.error
        if snps is not None and error.size: error = error[snps, :]
        sample_set = set(samples)
        sample_index_map = util.dict_invert(dict(enumerate(p.sample_index)))
        ibd = im.segment.SegmentSet(im.segment.Segment(x.snp, map(lambda y: (sample_index_map[y[0]], y[1]), x.samples),
                                                       x.bp, error_snps=x.error_snps) for x in self.info.ibd if (sample_set >= set([y[0] for y in x.samples])))
        info = ProblemInfo(p, g, snp=(self.info.snp[snps] if snps is not None else self.info.snp), ibd=ibd)
        return Problem(p, g, haplotype=h, info=info, error=error, frames=frames, lam=self.lam)
        
    def remove_nodes(self, nodes):
        '''A convenience method that returns a sub-problem without the nodes 'nodes'.'''
        return self.sub_problem(np.setdiff1d(self.pedigree.graph.nodes(), nodes))

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def info(self):
        '''Auxiliary phasing info struct.'''  
        return self._info

    @info.setter
    def info(self, info):
        '''Auxiliary phasing info struct.'''  
        self._info = info

    @property
    def num_errors(self):
        '''Genotype error counter.'''  
        return len(np.where(self.error != im.constants.OK)[0])

    @property
    def num_pedigree_samples(self):
        '''# people in pedigree.'''  
        return self.pedigree.n

    @property
    def num_data(self):
        '''Total # of data entries.'''  
        return self.genotype.num_data

    @property
    def num_genotypes(self):
        '''#genotype entries.'''  
        return self.genotype.num_samples * self.genotype.num_snps

    @property
    def num_samples(self):
        '''# genotyped people.'''  
        return self.genotype.num_samples

    @property
    def num_snps(self):
        '''#SNPs for each person.'''  
        return self.genotype.num_snps

    @property
    def snp_range(self):
        '''Find and cache the range of all SNPs.'''
        if self._snp_range is None:
            self._snp_range = np.arange(0, self.num_snps)
        return self._snp_range
    
    @property
    def components(self):
        '''Return the tuple (genotype, haplotype).'''  
        return self.genotype, self.haplotype
    
    @property
    def data(self):
        '''Return the tuple (genotype.data, haplotype.data).'''
        return self.genotype.data, self.haplotype.data

    @property
    def g(self):
        '''Return the genotype.data array.'''
        return self.genotype.data

    @property
    def h(self):
        '''Return the haplotype.data array.'''
        return self.haplotype.data

    @property
    def first_family(self):
        '''Return the first family in this object.'''
        return self.families(genotyped=False)[0]
    
####################################################################################
class ProblemInfo(object):
    '''A struct holding auxiliary information acquired during phasing: SNP annotation,
    recombination locations (in families for now), and a global IBD dictionary.  
    This is a lightweight stand-alone object that can be efficiently pickled and unpickled by
    partial phasing runs.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, pedigree, genotype, family_info=None, snp=None, ibd=None):
        '''Initialize from pedigree and genotype data sets.'''
        self.num_snps = genotype.num_snps
        self.num_samples = genotype.num_samples
        self.num_pedigree_samples = pedigree.n
        self.family_info = family_info if family_info else {}
        self.ibd = ibd if ibd else im.segment.SegmentSet() 
                    
        # SNP metadata            
        if snp is not None:
            # Import passed-in annotations
            self.snp = snp
        else:
            # Allocate record array containing SNP annotation        
            self.snp = np.zeros((self.num_snps,),
                             dtype=[
                                    ('chrom', np.uint8),  # Chromosome # containing the SNP
                                    ('name', np.chararray),  # SNP name (e.g., 'rs...')
                                    ('dist_cm', np.float),  # Genetic position
                                    ('base_pair', np.uint),  # Base pair position on chromosome
                                    ('count', '(4,)i4'),  # Genotype count: (1,1),(1,2),(2,2)
                                    ('frequency', '(4,)f4'),  # Genotype frequency: (1,1),(1,2),(2,2)
                                    ])
            # Import annotations from Genotype object
            for col in ('chrom', 'name', 'dist_cm', 'base_pair'):
                self.snp[col] = genotype.snp[col]

    #---------------------------------------------
    # Operators
    #---------------------------------------------
    def __key(self):
        '''Uniquely-identifying key of this object.'''
        return (self.num_snps, self.num_samples, self.num_pedigree_samples, self.family_info)

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
        return 'ProblemInfo[snps=%d, samples=%d, pedigree=%d]' \
            % (self.num_snps, self.num_samples, self.num_pedigree_samples) 
               
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def set_family_info(self, family, info):
        '''Family-to-info dictionary entry setting.'''
        if not self.family_info.has_key(family):
            self.family_info[family] = {im.constants.PATERNAL: None, im.constants.MATERNAL: None}
        self.family_info[family][info.parent_type] = FamilyInfo.from_info(info)
        
    def allele_frequency(self, allele):
        '''Estimated allele frequencies at each SNP (allele=1,2).'''
        count = self.snp['count']
        # If columns are 11,12,22 then they correspond to p^2,2pq,q^2 and p=p^2+0.5*2*p*q. Similarly q.
        return (2 * count[:, (0 if allele == 1 else 2)] + count[:, 1]) / (2.0 * np.sum(count[:, 0:3], axis=1))

    def snp_by_name(self, snp_name):
        '''Return the SNP index corresponding to a SNP name.'''
        return np.where(self.snp['name'] == snp_name)[0]
    
    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def genotype_frequency(self):
        '''Estimated genotype frequencies at each SNP (cols correspond to the genotypes
        (1,1),(1,2),(2,2)).'''  
        return self.snp['frequency']
    
    @property
    def snp_range(self):
        '''Return the range of all SNPs.'''
        return np.arange(0, self.num_snps)

    @property
    def pprint(self):
        return ProblemInfo.Printer(self.snp)
    
    class Printer(object):
        def __init__(self, snp):
            self.snp = snp
            
        def __getitem__(self, *args):
            '''Pretty-print SNP information.'''
            return ProblemInfo.Printer(self.snp.__getitem__(args[0])) if args else self

        def __repr__(self):
            s = ''
            for (a, b, c) in self.snp[np.array(['name', 'count', 'frequency'])]:
                alleles = 2.0 * sum(b[0:3])
                s += ('%-15s' % (a,) + ' ' + util.cprint(b, '%5d') + '  ' + util.cprint(c, '%7.5f')
                      + ' [%7.5f %7.5f]' % ((2 * b[0] + b[1]) / alleles, (2 * b[2] + b[1]) / alleles))
                s += '\n' 
            return s                

####################################################################################
class FamilyInfo(object):
    '''Holds auxiliary information acquired about a nuclear family: errors, recombination locations, etc.
    This is a lightweight stand-alone object that can be efficiently pickled and unpickled by
    partial phasing runs.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, template, recombinations):
        '''Initialize from fields. Note: recombinations is a SNP *original-ID-in-problem*
        (not index) array.'''
        # self.snps = snps
        self.template = template
        self.recombinations = recombinations

    @staticmethod
    def from_info(info):
        '''Convert a FamilyIbdComputer object into a FamilyInfo object.'''
        return FamilyInfo(info.template, info.recombination_snp)

    def copy(self):
        '''Deep copy.'''
        return FamilyInfo(self.template, self.recombinations.copy())

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __key(self):
        '''Uniquely-identifying key of this object.'''
        return (self.template, self.recombinations.tolist())

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
        return 'FamilyInfo[template=' + repr(self.template) + ']'

####################################################################################
class LambdaEvaluator(object):
    '''A functor that evaluates an interpolant lambda(f) given discrete estimates from the samples.''' 
    def __init__(self, data):
        '''Load data array with the columns f, lam.'''
        self.__F = data[:, 0]
        self.__L = data[:, 1]

    def __call__(self, fvals):
        '''Return an interpolant lambda(f) at fvals given the discrete binned values F, L obtained by
        lambda_mean().'''
        return np.interp(fvals, self.__F, self.__L)
