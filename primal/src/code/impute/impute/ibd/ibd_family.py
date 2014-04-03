'''
============================================================
IBD segment identification between a proband and its
gentoyped parent based on recombination event estimation
between the relevant haplotypes.

This is a self-sufficient implementation of the IBD segment
abstract factory.

Created on July 22, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
from impute.data.constants import INDETERMINATE, PATERNAL, ALLELES
from impute.tools import genotype_tools as gt
from impute.ibd import ibd
from impute.tools.param import PhaseParam

#---------------------------------------------
# Methods
#---------------------------------------------
def ibd_segments_parent_child(haplotype, parent, child, parent_type, het_snps,
                              min_segment_length=INDETERMINATE, debug=False):
    '''Return 1) IBD segments between a parent and child. 2) List of het_snp indices at which there
    are likely genotype errors.'''
    # Parent phase is dummy. Could be either one.
    # if debug:
#    if parent == 1036 or child == 1036:
#        print 'parent-child IBD', parent, child, parent_type
    return ibd.ibs_segments(haplotype, parent, child, PATERNAL, parent_type, het_snps, debug=debug)

def confidence(snp): raise ValueError('To be implemented')

def ibd_segments_in_family(h, family, parent_het_fill_threshold, debug=False):
    '''A generator helper method that yields IBD segments in a family (compare each parent
    with the corresponding child haplotype and output segments where the haps are IBS).
    h=haplotype data set. family=family to process. parent_fill_threshold=% of parent haplotypes
    required to be phased to use it to find segments.'''
    for parent_type, parent in family.parents_dict().iteritems():
        # print family, parent, h.fill_fraction(sample=parent)
        if h.fill_fraction(sample=parent) > parent_het_fill_threshold:
            het_snps = gt.where_heterozygous(h.data, parent)
            for child in family.children:
                yield ibd_segments_parent_child(h, parent, child, parent_type, het_snps, debug=debug)

def ibd_segments_in_duo(h, parent, child, parent_type, parent_het_fill_threshold, debug=False):
    '''Similar to ibd_segments_in_family, but for a parent-child duo.'''
    if h.fill_fraction(sample=parent) > parent_het_fill_threshold:
        het_snps = gt.where_heterozygous(h.data, parent)
        yield ibd_segments_parent_child(h, parent, child, parent_type, het_snps, debug=debug)

def ibd_segments_in_family_sibs(request, family):
    '''Even if both parents are not genotyped, output IBD segments within siblings, since
    they share (long) segments on their corresponding paternal haplotypes and maternal haplotypes.
    Only long segments are considered.'''
    problem, params = request.problem, request.params
    genotyped_children = gt.genotyped_children(problem, family)

    # TODO: swap the haps of phased children (that are parents in another family) if inconsistent with their sibs       
    # self.d = (self.h_child != np.transpose(np.tile(self.h_child[:,template], (self.num_children,1)))).astype('int8')
        
    # There exist at least two genotyped children to compare. Pick one as a template
    # and compare every other child against it, or all children as templates
    h = problem.haplotype
    template_mode = params.template_mode
    if template_mode == PhaseParam.TEMPLATE_MODE.ONE:
        # 1 templates = most-filled child 
        fill = problem.fill_fraction(sample=genotyped_children)
        template_index = np.argmax(fill[:, 1])
        template_child = int(fill[template_index, 0])
        templates = [(template_index, template_child)]
        pairs_to_check = lambda i, j: i != j
    if template_mode == PhaseParam.TEMPLATE_MODE.TWO:
        # 2 templates = two most-filled children 
        fill = problem.fill_fraction(sample=genotyped_children)
        template_index = np.argsort(-fill[:, 1])[0:2]
        template_child = fill[template_index, 0].astype(np.int8)
        templates = [(x, template_child[i]) for (i, x) in enumerate(template_index)]
        pairs_to_check = lambda i, j: i != j
    elif template_mode == PhaseParam.TEMPLATE_MODE.ALL:
        # EVery child serves as a template
        templates = enumerate(genotyped_children)
        pairs_to_check = lambda i, j: i < j
    else:
        raise ValueError('Unsupported template mode ''%s''' % (template_mode,))
        
    for parent_type in ALLELES:
        for i, template_child in templates:
            for j, child in enumerate(genotyped_children):
                if pairs_to_check(i, j):
                    yield ibd.ibs_segments(h, template_child, child, parent_type, parent_type,
                                           error_filter='fast_max',
                                           min_segment_length=params.min_segment_length_sibs,
                                           include_alt_phase=False, debug=params.debug)

#def align_gender_phases(request, family):
#    '''If both parents are not genotyped, we need to align the gender phases of the children.
#    That is, make sure that hap 0 is paternal and hap 1 is maternal (or vice versa) in all children.
#    Algorithm:
#    - Find all IBD segments using GERMLINE.
#    - Find '''
#    problem, params = request.problem, request.params
#    genotyped_children = gt.genotyped_children(problem, family)
    
#---------------------------------------------
# Private Methods
#---------------------------------------------
