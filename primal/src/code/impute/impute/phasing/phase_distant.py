'''
============================================================
Phase using IBD segment identification between distant
individuals.

Created on August 22, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, impute as im, itertools as it
from chain import Filter
from impute.ibd import ibd
from impute.phasing.phase_core import new_phaser_chain, FamilyPhaser
from impute.ibd.distant.ibd_distant import RelativeCollection
from impute.ibd.distant.ibd_hmm import prob_ibd_hmm
from impute.ibd.segment import START, STOP
from impute.phasing.pre_processing import phased_samples_phaser
from impute.data.constants import PATERNAL, MATERNAL

#---------------------------------------------
# Methods
#---------------------------------------------
def phase_by_segment_priority(problem, segments):
    '''Given a list of segments, determine the best segment to use at each SNP, and phase using
    the max consensus. This assumes that segments are between pairs of haplotypes, one of which
    is filled.'''
    problem.info.ibd += segments  # Not grouping to disjoint here; TODO: maybe in the future?
    # Prioritize which segment to use for phasing at each SNP
    s = im.idist.best_segment_index(problem, segments)
    # For each segment, use it to phase at all SNPs at which it is the best segment
    for i, segment in enumerate(segments):
        ibd._phase_in_ibd_segment(problem, np.where(s == i)[0], np.array(list(segment.samples), dtype=np.uint), 'max')

####################################################################################
def __handle_sib_comparison(self, request, family):
    '''Locate child-child IBD segments in families none of the parents are genotyped.'''
    # h = problem.haplotype
    problem, params = request.problem, request.params
    genotyped_children = [x for x in family.children_list if problem.is_genotyped(x)]
    if not np.any([problem.is_genotyped(parent) for parent in family.parents]) and len(genotyped_children) >= 1:
        if params.debug:
            print '=' * 80
            print 'Running __handle_sibs() in family', family, 'genotyped_children', genotyped_children
            print '=' * 80
        __handle_sib_founder_family(self, request, family)
    
def __handle_sib_founder_family(self, request, family):
    '''A) Align the POO phases of all phased quasi-founder siblings. 
    
    B) If alignment was successful, phase unphased sibs in the same family using IBD segments
    between them and the phased siblings.
     
    Output long IBD segments among the siblings. We only output long segments since we know
    sibs share long segments with the parents.'''
    problem, params = request.problem, request.params
    g, h = problem.g, problem.haplotype
    genotyped_children = im.gt.genotyped_children(problem, family)
    
    # Find phased children, identify IBD segments, build paternal haplotypes (by coloring
    # child haplotypes), align children POO phases
    phased_children = np.array([x for x in genotyped_children if h.fill_fraction(sample=x) >= params.surrogate_parent_fill_threshold])
    if params.debug:
        print 'Phased children', phased_children
    if not phased_children.size:
        # No phased children
        return
    elif phased_children.size == 1:
        # If there's just one kid, no need to align
        poo, separation = np.array([1.0]), 1.0
    else:
        # At least two children ==> look for segments and try to align their POO phase.
        segments_phased = im.ibd_distant_hap.among_samples_segments(problem, phased_children, request.params)
        if params.debug:
            print 'Segments among phased children:'
            print segments_phased
        if segments_phased.length:
            pa = im.color.hap_color.hap_colors(list(it.product(sorted(phased_children), im.constants.ALLELES)), segments_phased, max_colors=4)
            poo, separation, _, _ = im.color.hap_color.best_hap_alignment_to_colors(pa)
            if params.debug:
                print pa
                print 'POO phases', poo
                print 'Separation measure', separation 
                print 'Parental haplotype coverage', pa.color_sequence_coverage(np.arange(4))
        else: separation = 0.0 # no segments, can't align and can't trust kids phased the unphased samples 
            
    if np.abs(separation) > params.poo_coloring_measure_threshold:
        # Alignment of all kids succeeded
        if params.debug:
            print 'Alignment successful'
        
        # Align POO phases - flip haplotypes of children with reversed haplotype
        # This is a -local- alignment within this family. Still need to globally align all POO phases.
        # This is done by the poo module, after phasing and IBD segment calculation+index are done.
        flipped_children = phased_children[poo < -params.poo_coloring_measure_threshold]
        h_children = problem.h[:, flipped_children, :]
        problem.h[:, flipped_children, PATERNAL] = h_children[:, :, MATERNAL]
        problem.h[:, flipped_children, MATERNAL] = h_children[:, :, PATERNAL]
        
        # Phase each unphased sib 'sample' using IBD segments between sample and phase_children  
        # Note: phased sibs these are all full sibs of the proband, not half sibs. So they are
        # m=2 meioses away from the proband.
        for sample in [x for x in genotyped_children if h.fill_fraction(sample=x, snps=im.gt.where_heterozygous(g, x)) < params.het_fill_threshold]:
            relatives = RelativeCollection.from_sibs(sample, phased_children)
            if relatives.length > 0:
#                segments = sum((im.ibd_distant_hap.hap_segments(IbdProblem(problem, (sample, allele), (sib, allele), None, params)) 
                                # for sib, allele in it.product(phased_children, ALLELES)), im.segment.SegmentSet([]))
                segments = im.idist.ibd_segments_with_relatives(problem, sample,
                                                                relatives.info['index'], params,
                                                                prob_ibd_calculator=prob_ibd_hmm)
                # Since j is POO-aligned and is a sib of i, segments can only be between (i,a),(j,a).
                # Replace the dummy i allele in segments by j's allele. 
                segments = im.segment.SegmentSet(im.segment.Segment(x.snp, [(x.samples[0][0], x.samples[1][1]), x.samples[1]], x.bp,
                                                                    error_snps=x.error_snps, confidence=x.confidence, cm=x.cm)
                                                 for x in segments)
                if params.debug:
                    print 'Segments after assigning sib alleles to proband haplotypes:'
                    print segments
                phase_by_segment_priority(problem, segments)
    else:
        if params.debug:
            print 'Alignment unsuccessful'


####################################################################################
class __ParentDistantPhaser(Filter):
    '''Determine IBD segments between an unphased individual and its distant (but not too distant)
    relatives using the individual's trivially-phased homozygous SNPs.'''

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __init__(self, next_filter=None, name=None, single_sample=None, max_samples=np.inf, debug=False):
        Filter.__init__(self, name=name, next_filter=next_filter)
        # Debugging options
        self.single_sample = single_sample
        self.debug = debug
        self.max_samples = max_samples
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def handle(self, request):
        '''Handle the phasing request: perform multiple pass. In each pass, loop over unphased individuals and
        compare each to a set of its distant relatives on the undirected pedigree graph.
        Use HMM to estimate IBD segments.'''
        for phased_fill, max_path_length, target_fill in request.params.distant_phasing_params:
            self.__pass(request, phased_fill, max_path_length, target_fill)

    def __pass(self, request, phased_fill, max_path_length, target_fill):
        '''A single pass of distant relative phasing : loop over unphased individuals and
        comSpare each to a set of its distant relatives on the undirected pedigree graph.
        Use HMM to estimate IBD segments.'''
        problem, params = request.problem, request.params
        # debug = params.debug
        debug = True
        chunk_size = 2  # Chunk size of relatives to process at a time
        chunk_growth_size = 1.5  # 1.2  # Chunk size growth factor
        if debug:
            print '/' * 80
            print 'Distant phasing pass: phased_fill %.2f, max_path_length %d, target_fill %.2f' % \
            (phased_fill, max_path_length, target_fill)
            print '/' * 80 
        if self.single_sample is not None:
            samples = [self.single_sample] 
            filled = problem.fill_fraction()
        else:
            filled = problem.fill_fraction()
            samples = filled[np.where(filled[:, 1] < phased_fill)[0], 0].astype(int)
            filled = filled[samples, :]
            filled = filled[np.argsort(filled[:, 1])]
            samples = filled[:, 0].astype(int)
            if params.selected_mode:
                # Restrict samples to selected samples only. Use a boolean array
                # 'selected' to indicate which entries of 'samples' are in the selected set.
                selected = np.in1d(params.selected_samples, samples)
                samples, filled = samples[selected], filled[selected]
        if debug:
            print 'Target phasing fraction: %.2f, max_path_length %d' % (target_fill, max_path_length)
            print 'samples to be phased (< %.2f phased):' % (phased_fill,) 
            print filled
        for i, sample in enumerate([int(x) for x in samples]):            
            if debug:
                print '#' * 60
                print 'Distant phasing sample %d' % (sample,)
                print '#' * 60
            # Find all relatives and sort them by sc depth, then by desc fill%
            relatives = RelativeCollection.in_neighborhood(sample, problem, 1, max_path_length, phased_fill)
            if relatives.length:
                fill = np.concatenate((problem.fill_fraction(sample=relatives.info['index']), relatives.info['distance'][:, 0][np.newaxis].transpose()), axis=1)
                index = np.lexsort((-fill[:, 1], fill[:, 2]))
                sorted_relatives = relatives.info['index'][index]
                # Define a sequence of increasingly-larger chunks. TODO: wrap in iterator
                n = relatives.length
                ind = np.concatenate(([0], np.round(chunk_size * chunk_growth_size ** 
                                                    np.arange(np.floor(np.log(n) / np.log(chunk_growth_size) 
                                                                       - np.log(chunk_size)))).astype(int)))
                if ind[-1] < n: ind = np.concatenate((ind, [n]))
                if debug: print 'Chunk size ind', ind
            else:
                sorted_relatives = []
                ind = []
            
            # Look for increasingly-distant relatives until we get IBD coverage sufficient to phase
            # the entire sample's chromosome 
            fill = problem.fill_fraction_of_sample(sample)
            k = 0
            for k, chunk in enumerate(map(lambda k: (ind[k], ind[k + 1]), range(len(ind) - 1))):
                relatives_chunk = sorted_relatives[chunk[START]:chunk[STOP]]
                if debug:
                    print '-' * 70
                    print 'Processing sample %d, %d/%d fill %f, chunk #%d' % \
                    (sample, i + 1, len(samples), fill, k + 1)
                    print 'Relatives fill'
                    print problem.fill_fraction(sample=relatives_chunk)
                    if self.single_sample is not None: print 'h[0,%d] = %s' % (self.single_sample, repr(problem.h[0, self.single_sample, :]))
                if relatives.length:
                    segments = im.idist.ibd_segments_with_relatives(problem, sample, relatives_chunk, params,
                                                                    prob_ibd_calculator=prob_ibd_hmm)
                    # Not grouping to disjoint at this point, maybe in the future
                    if debug:
                        print 'IBD segments, sample-surrogate parent'
                        print segments
                        print 'Phasing sample by segment priority'
                    phase_by_segment_priority(problem, segments)
                # Increment relative set depth
                fill = problem.fill_fraction_of_sample(sample)
                if fill >= target_fill:
                    if debug: print 'Reached target at sample %d, %d/%d fill %f' % (sample, k + 1, len(samples), fill)
                    break
            if debug: print 'Done with sample %d, %d/%d fill %f' % (sample, k + 1, len(samples), fill)
            # Debugging, process only the first max_samples samples
            if i == self.max_samples - 1: break
        return False

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __phase_distant(self, problem, segments, sample):
        s = im.idist.best_segment_index(problem.info.snp['base_pair'], segments)
        for i, segment in enumerate(segments):
            ibd._phase_in_ibd_segment(problem, np.where(s == i)[0], np.array(list(segment.samples)), 'max')

####################################################################################
'''Main phasing processing -- distant IBD.'''
def family_sib_comparison_phaser(next_filter=None, debug=False, print_times=False):
    chain = [
             FamilyPhaser(family_handle=__handle_sib_comparison, name='* Sibs', genotyped_families=False),
             phased_samples_phaser  # In selected samples mode, reset all non-selected haplotypes back to the original values 
             ]     
    return new_phaser_chain(chain, name='Family Sib', next_filter=next_filter,
                            debug=debug, print_times=print_times)

def distant_phaser(next_filter=None, debug=False, print_times=False,
                   max_samples=np.inf, single_sample=None):
    chain = [
             __ParentDistantPhaser(name='* Distant IBD',
                                   max_samples=max_samples,
                                   single_sample=single_sample,
                                   debug=debug),
             phased_samples_phaser  # In selected samples mode, reset all non-selected haplotypes back to the original values 
             ]    
    return new_phaser_chain(chain, name='Distant IBD', next_filter=next_filter,
                            debug=debug, print_times=print_times)
