'''
============================================================
Tools for estimating IBD segments between children in
a nuclear family.

Created on July 6, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import util, numpy as np, impute as im
from scipy import ndimage
from impute.data.constants import PATERNAL
from impute.ibd.segment import SegmentSet, Segment
from impute.ibd import ibd, segment, diff
from impute.tools import genotype_tools as gt

#---------------------------------------------
# Methods
#---------------------------------------------
class ChildComparator:
    '''Compares children haplotypes and one parent in a nuclear family. Outputs recombination
    events as a FamilyIbdComputer object.'''
    
    def __init__(self, request, family, children=None):
        '''Initialize an IBD segment computer for a particular family in a problem. A children set
        of interest must be specified (typically, those with high fill %).'''
        self.request = request
        self.problem = request.problem
        self.params = request.params
        self.family = family
        # Children to consider
        self.children = family.children_array if children is None else children
        if self.params.debug:
            print 'Filled children', self.children
        # self.filled_children = children         

    def phase_parent_by_template(self, info):
        '''Impute the parent haplotype. A random parent random phase is selected and switched at
        each recombination location in the template child.'''
        bp = self.problem.info.snp['base_pair']
        num_snps = self.problem.num_snps
        parent = info.parent
        child = info.template_child
        parent_type = info.parent_type
        if self.params.debug:
            print 'phase_parent_by_template(parent=%d, template=%d, parent_type=%d)' % (parent, child, parent_type)
        # h = problem.haplotype.data                 
        # Infer the parent's haplotypes from the template child's recombinations r
        r = info.recombination_snp
        edge = [y for (x, y) in r if x == child]
        # r = np.concatenate(([-1], edge, [problem.num_snps-1]))
        
        phase = PATERNAL  # Random phase set at the starting of the chromosome
        segments = segment.edges_to_segments(self.problem.snp_range, edge, phase,
                                             self.problem.num_snps, cover=True)
        segment_set = SegmentSet((Segment((x[0], x[1]), ((parent, x[2]), (child, parent_type)),
                                          (bp[x[0]], im.segment.stop_bp(bp, x[1], num_snps))) 
                                  for x in segments))
        # Only a single sweep suffices here, since only two samples are involved. If a homogeneous
        # entry is discovered and used to propagate phases in the second of the two IBD sets in
        # a SNP range, there's no need to fix the first one since the data is already fully filled.
        ibd.phase_by_ibd(self.request, [segment_set], 'max', num_sweeps=1)
            
    def child_recombinations(self, parent_type, template_child=None, remove_errors=True):
        '''A helper function to handle_founder_parent_by_children() that processes a single family.
        Returns the list of recombinations in the children of a family w.r.t. to the haplotypes of
        a single parent (parent_type=PATERNAL/MATERNAL). template_child the child to be used as a template
        against which all other children are compared. If not specified, the first child on the list is used.
        
        This method will modify the problem object by flagging and removing errors if and only if
        remove_errors = True.'''
        # Initializate auxiliary object holding recombination
        info = FamilyIbdComputer(self.problem, self.params, self.family, parent_type, self.children)
        if not info.snps_exist:
            return ([], [], info)
        
        template_index = np.where(info.children == template_child)[0][0] if template_child is not None else 0
        info.set_template(template_index)
        
        # Find genotype errors
        errors = info.find_errors()
        # Flag genotype errors in the problem object. TODO: move to the standard phase_by_ibd() call?
        if remove_errors:
            for (error_type, (snp, child)) in errors.iteritems():
                msg = FamilyIbdComputer.ERROR_MESSAGES[error_type]
                self.problem.genotype_error(snp, child, 'Children comparison spike - ' + msg)
        
        # Find recombinations
        info.find_recombinations()
        # info.plot()
        # plt.savefig('family_filtered.png')
        # print 'errors', len(errors), errors
        return ([], errors, info)
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------

####################################################################################
class FamilyIbdComputer(object):
    '''A helper struct containing family information to pass around the unphased-parent-from-
    children phasing internal loop.'''
    INDECISIVE = -1
    TEMPLATE = 0
    NON_TEMPLATE = 1
    
    ERROR_MESSAGES = {NON_TEMPLATE: 'non-template children',
                   TEMPLATE: 'template child',
                   INDECISIVE: 'indecisive agreement'}
    ERROR_TYPES = ERROR_MESSAGES.keys()
    
    def __init__(self, problem, params, family, parent_type, children, debug=False):
        self.problem = problem
        self.params = params
        self.sample_id = self.problem.pedigree.sample_id
        self.sample_index = self.problem.pedigree.sample_index
        (self.g, self.h) = self.problem.data
        self.family = family
        self.parent_type = parent_type
        self.parent = family.parents[parent_type]
        # self.children = np.array(list([x for x in family.children if problem.is_genotyped(x)]))#[np.array([0,1,2,5,6,7,8,9,10,11,12])]
        self.children = children
        self.num_children = len(self.children)
        self.num_snps = self.problem.num_snps
        self.snps = None
        self.het_snps = gt.where_heterozygous(self.g, sample=self.parent)
        self.template = -1
        self.__filter_snps()
        self._recombination_index = np.empty((0, 2), dtype=int)
        self.error_filter_length = 5
        self.debug = debug
    
    def __repr__(self):
        return '[children=%s  parent=%d  type=%d  #relevant snps=%d  template=%d]' % \
            (repr(self.children), self.parent, self.parent_type, len(self.snps), self.template)
        
    def __filter_snps(self):
        '''Find SNPs that are heterozygous in the parent and filled in all children haplotypes.
        Sets the h_child (child haplotype array) and snps fields of this object.'''
        # Find het SNPs in parent
        
        # Further filter SNPs so that all children have data
        h_child = self.h[:, self.children, self.parent_type][self.het_snps, :]
        filled = np.where(np.min(h_child, axis=1) > 0)[0]
        # filled = np.where(np.sum(h_child > 0, axis=1) >= 5)[0]
        self.snps_exist = filled.size > 0
        if self.snps_exist:
            self.h_child = h_child[filled, :]
            self.snps = self.het_snps[filled]

    def set_template(self, template):
        '''Compute differences d between template child and all others.'''
        self.template = template
        self.d = (self.h_child != np.transpose(np.tile(self.h_child[:, template],
                                                       (self.num_children, 1)))).astype('int8')
        
    def find_errors(self):
        '''Detect spikes in d = likely genotype errors. A spike in a single child
        means an error in that child. A spike in all but the template means
        an error in the template. The other indecisive cases are marked as errors.''' 
        
        # Use median filter to detect spikes in d
        filtered_diff = ndimage.median_filter(self.d, (self.error_filter_length, 1))
        # filtered_diff = ndimage.median_filter(self.d, (self.error_filter_length,1), mode='constant', cval=0)
        # ndimage.median_filter(self.d, (3,1), mode='constant', cval=0)
        # filtered_diff = util.max_except_center_filter1d(self.d, 5) 
        # print self.d
        # print filtered_diff
        # print filtered_diff2
        # filtered_diff = util.max_except_center_filter(d_snps, error_filter_length)
        difference = (self.d != filtered_diff).astype('int8')
        num_diff = np.sum(difference, axis=1)
        errors = np.nonzero(num_diff)[0]
        error_dict = {}
        if errors.size:
            num_errors = num_diff[errors] if errors.size else np.array([])
            for error_type in FamilyIbdComputer.ERROR_TYPES:
                (snp_index, child_index) = self.__error_index(difference, errors, num_errors, error_type)
                error_dict[error_type] = (snp_index, child_index)

        # Remove errors from data array
        self.dorig = self.d.copy()
        self.d = filtered_diff
        return error_dict
        
    def find_recombinations(self):
        '''Set r._recombination_index to an array [child, event] of recombination events.
        event denotes an edge index in the d array.
        - A change in a single child indicates recombination in that child.
        - A change in all children except the template implies a recombination in the template. 
        - All other cases are indecisive.''' 
        edge = np.diff(self.d, axis=0)
        # deriv = ndimage.convolve1d(self.d, [1, 0, -1], axis=0)
        
        # Number of children an edge occurred in 
        num_edges = np.sum(edge != 0, axis=1)
        
        r = []
        # Non-template children edges
        location = np.where(num_edges == 1)[0]
        if location.size:
            child_index = self.children[np.nonzero(edge[location, :])[1]]
            r += zip(child_index, location)
        
        # Template child edges
        location = np.where(num_edges == self.num_children - 1)[0]
        if location.size:
            child_index = np.tile(self.children[self.template], (len(location),))
            r += zip(child_index, location)

        # Indecisive cases: shouldn't happen!
        indecisive = np.where(np.logical_and(num_edges != 0,
                                             np.logical_and(num_edges != 1,
                                                            num_edges != self.num_children - 1)))[0]
        if indecisive.size:
            # Indecisive cases encountered
            # log(WARN, 'Indecisive recombination situations encountered. Did you filter genotype errors first?')
            # log(WARN, 'indecisive ' + repr(self.snps[indecisive]) + ' num_edges ' + repr(num_edges[indecisive]))
            # raise ValueError('Indecisive recombination situations encountered. Did you filter genotype errors first with find_errors()?')
            location = np.nonzero(edge[indecisive, :])
            r += zip(self.children[location[1]], indecisive[location[0]])
        r = np.array(r)
        # Sort recombinations by child ID, then by SNP location
        self._recombination_index = r[np.lexsort((r[:, 1], r[:, 0]))] if r.size else r

    def parent_fill_fraction(self):
        return self.problem.haplotype.fill_fraction(sample=self.parent)
        
    def ibs_segments(self):
        '''Convert recombinations r to IBD segments. A generator.'''
        if not self.snps_exist:
            return
        r = self.recombination_index
        bp = self.problem.info.snp['base_pair']
        # print r
        
        # Calculate IBD segments between all children and parent 
        for child in self.children:
            # print '-------------------- child %d ----------------------' % (child,)
            edge = [y for (x, y) in r if x == child]
            phase = PATERNAL  # Random phase set at the starting of the chromosome
            segments = segment.edges_to_segments(self.snps, edge, phase, self.num_snps, cover=True)
            if self.debug:
                print 'child', child, 'edge', edge, 'segments', segments

            # Flip initial phase if parent is phased at hap snps by comparing its 
            # haplotypes with the child haplotypes on the longest IBD segment
            i = np.argmax(np.diff(segments)[:, 0])
            segment_phase = phase if np.mod(i, 2) == 0 else 1 - phase
            actual_phase = self.__equal_parent_hap(child, segments[i, 0], segments[i, 1])
            if segment_phase != actual_phase:
                segments = segment.flip_phase(segments)

            # Output segments in the standard IBD segment format
            yield SegmentSet((Segment((x[0], x[1]), ((self.parent, x[2]), (child, self.parent_type)),
                                      (bp[x[0]], im.segment.stop_bp(bp, x[1], self.num_snps)))
                              for x in segments))

#    def plot(self, title=None, template=True, y_labels='index'):
#        '''Create haplotype coloring plot.'''
#        if y_labels == 'index':
#            y = self.sample_index[self.children]
#        elif y_labels == 'id':
#            y = self.sample_id[self.children]
#        else:
#            y = self.children 
#        t = self.template_child if template else self.parent
#        index = self.children if template else np.concatenate((np.array([self.parent]), self.children))
#        return plots.plot_comparison(template=self.h[:,t,self.parent_type],
#                                     haps=self.h[:,index,self.parent_type],
#                                     snps=self.snps,
#                                     title=title, y=y)
        
    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def template_child(self):
        '''Return the template child ID.'''
        return self.children[self.template]
    
    @property
    def recombination_index(self):
        '''Return recombination indices in the snps array.'''
        return self._recombination_index
    
    @property
    def recombination_snp(self):
        '''Return recombination SNP locations.'''
        r = self._recombination_index.copy()
        if r.size:
            r[:, 1] = self.snps[r[:, 1]]
            return r
        else:
            return r
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __equal_parent_hap(self, child, snp_start, snp_stop):
        '''Return the list of equal haplotypes for a snp_segment=(snp_start,snp_stop)
        in the samples id1,id2.
         
        Haplotype comparison is delegated to self.hap_comparator and must be at least
        self.threshold for a pair to be considered equal.
        
        The function emits tuples (i,j), where i and j are allele indices for which
        hapotypes id1-i and id2-j are equal.'''        
        h = self.h
        index = np.where([diff.hap_comparator_difference(h[snp_start:snp_stop + 1, self.parent, j],
                                                         h[snp_start:snp_stop + 1, child, self.parent_type]) 
                          >= self.params.prob_ibd_threshold for j in (0, 1)])[0]
        return index[0] if index.size else None 

    def __error_index(self, diff, errors, num_errors, error_type):
        '''Return the corresponding snp and child indices of genotype errors. These are rows
        in the errors array that have error_type non-zeros.'''
        if error_type == FamilyIbdComputer.NON_TEMPLATE:
            # Non-template children errors
            snp_index = errors[np.where(num_errors == 1)[0]]
            return (self.snps[snp_index], self.children[np.where(diff[snp_index, :] != 0)[1]])
        if error_type == FamilyIbdComputer.TEMPLATE:
            # Template child errors
            snp_index = errors[np.where(num_errors == self.num_children - 1)[0]]
            return (self.snps[snp_index], np.tile(self.children[self.template], (len(snp_index),))) 
        else:
            # Indecisive cases: flag parent+all children as errors -- the best we can do for now.
            snp_index = errors[np.where(np.logical_and(num_errors != 1, num_errors != self.num_children - 1))[0]]
            a = util.flattened_meshgrid(np.concatenate((self.children, [self.parent])), self.snps[snp_index])
            return (a[1], a[0])
