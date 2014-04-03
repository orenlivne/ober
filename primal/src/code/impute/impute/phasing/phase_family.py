'''
============================================================
Phase nuclear families by comparing parent and children
haplotypes to determine recombination points, collocating
parent and children haplotypes into IBD sets within segments,
and collectively using the elements of each such set to
phase all participating individuals.

Created on July 6, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import util, numpy as np
from impute.tools import pedigree_tools as pt 
from impute.ibd import ibd_family as ip
from impute.ibd import ibd_child as ic
from impute.ibd import ibd
from impute.tools import genotype_tools as gt
from impute.phasing.phase_core import FamilyPhaser, new_phaser_chain
from impute.data import constants
from chain import Filter
from impute.phasing.pre_processing import phased_samples_phaser

####################################################################################
def __handle_parent_child(self, request, family):
    '''Determine IBD segments between in a parent and a child => phase them based on
    each other (if one hap filled, use it to fill the other if missing).'''
    problem, params = request.problem, request.params
    return ibd.phase_by_ibd(request,
                            ip.ibd_segments_in_family(problem.haplotype, family,
                                                      parent_het_fill_threshold=params.het_fill_threshold,
                                                      debug=params.debug), 'max')

####################################################################################
def __handle_all_duos(self, request):
    '''Similar to __FamilyParentChildPhaser, only works with parent-child duos and does not compare
    multiple children in a family. Handles all duos; processes children in pre-ordering (top-to-bottom
    in the pedigree).'''    
    problem, params = request.problem, request.params
    # Build all duos 
    duos = np.empty((0, 3), dtype=int)
    for parent_type in constants.ALLELES:
        d = problem.duos(parent_type)
        duos = np.concatenate((duos, np.concatenate((d, np.tile(parent_type, (d.shape[0], 1))), axis=1)))
    # Pre-order nodes
    duos = duos[np.argsort([problem.pedigree.depth[x] for x in duos[:, 0]]), :]
    
    h = problem.haplotype
    for child, parent, parent_type in duos:
        ibd.phase_by_ibd(problem,
                         ip.ibd_segments_in_duo(h, parent, child, parent_type,
                                                parent_het_fill_threshold=params.het_fill_threshold,
                                                debug=params.debug),
                                 'max', debug=params.debug)
    return False

####################################################################################
def __handle_outer_duos(self, request):
    '''Similar to __FamilyParentChildPhaser, only works with parent-child duos and does not compare
    multiple children in a family. Handles only children outside nuclear families.'''
    problem, params = request.problem, request.params
    h = problem.haplotype
    family_members = problem.families_union(min_children=0)
    for parent_type in constants.ALLELES:
        for child, parent in pt.selected_duos(problem, params, parent_type):
            if not util.is_member(family_members, [child]):
                ibd.phase_by_ibd(request, ip.ibd_segments_in_duo(h, parent, child, parent_type,
                                                                 parent_het_fill_threshold=params.het_fill_threshold,
                                                                 debug=params.debug), 'max')
    return False

####################################################################################
def __handle_child_comparison(self, request):
    '''In families with at least 3 children:
    if one parent is a founder (or more generally, not sufficiently phased in het snps)
    and the other is not, the children in this family will be phased well by ParentChildFounder,
    but the non-founder parent will not be.
    
    By comparing children's haplotypes against a template child (=the most-filled child)
    and translating that into comparison between children haps and the unphased parent's, we can infer
    their IBS segments and subsequently the parent's haplotypes. 
    
    Note that the parent will have random hap-gender-assignment: we can't know which one of his/her
    haplotypes is paternal and which one is maternal (we might at a later stage, if his/her parent
    genotypes are genotyped or imputed by ancestor imputation).'''
    problem, params = request.problem, request.params 
    g, h = problem.components
    # Find families with at least min_consensus_samples genotyped children, or use single
    # family if debug mode (single_member) is on
    potential_families = problem.find_families_by_member(params.single_member, genotyped=False,
                                                         min_children=params.min_consensus_samples) \
    if params.single_member else pt.selected_families(problem, params, genotyped=False, min_children=params.min_consensus_samples)
    families = [f for f in potential_families
                if len(problem.find_samples_with_fill_ge(params.surrogate_parent_fill_threshold, sample=f.children_array))
                >= params.min_consensus_samples]
    if params.debug: print '__handle_child_comparison(), families to process', list(families)
    for family in families:
        genotyped_parent_dict = [(k, v) for (k, v) in family.parents_dict().iteritems() if problem.is_genotyped(v)]
        num_genotyped_parents = len(genotyped_parent_dict)
        # If both parents are genotyped, use all children - it is probably safe enough to generate
        # enough SNPs to work with (het in parent + filled in all children), since it has worked in the past.
        # If not both parents are genotyped, use filled children only to generate enough relevant SNPs.
        genotyped_children = np.array([x for x in family.children_array if problem.is_genotyped(x)])
        filled_children = genotyped_children if num_genotyped_parents == 2 else \
            problem.find_samples_with_fill_ge(params.surrogate_parent_fill_threshold, sample=genotyped_children)[:, 0].astype(np.int)
        comparator = ic.ChildComparator(request, family, filled_children)
        # for parent_type, parent in reversed(family.parents_dict().items()):
        for parent_type, parent in genotyped_parent_dict:
            # het_snps = gt.where_heterozygous(h.data, parent)
            het_snps = gt.where_heterozygous(g.data, parent)
            if h.fill_fraction(sample=parent, snps=het_snps) < params.het_fill_threshold:
            # if is_founder[parent]:
                # Choose template = most-filled child
                fill = problem.fill_fraction(sample=filled_children)
                if params.debug:
                    print '=' * 105
                    print 'Children comparison in', family, 'parent_type', parent_type
                    print '=' * 105
                    print [problem.is_genotyped(x) for x in family.children]
                    print '# genotyped children', sum(problem.is_genotyped(x) for x in family.children),
                    print '# parent het snps', len(het_snps)
                    print 'Filled children', filled_children
                    print 'Family''s fill:\n', problem.fill_fraction(sample=family.member_set)
                template_child = int(fill[np.argmax(fill[:, 1]), 0])
                if params.debug:
                    'template_child', template_child
                # Choose a template child at random (first index in the children list)
                # template_child = list(family.children)[0]
                (_, _, info) = comparator.child_recombinations(parent_type, template_child)
                # Save selected entries from the family info class in problem
                problem.set_family_info(family, info)
                # Impute parent from the template child
                if params.debug:
                    print family, parent_type
                    print 'Child recombinations'
                    print info.recombination_snp
                comparator.phase_parent_by_template(info)
                # Now phase children (and possibly some more of the parent) using IBD segments
                # found among them and the parent
                ibd.phase_by_ibd(request, info.ibs_segments(), 'majority')
    return False

#---------------------------------------------
# Private Methods
#---------------------------------------------

'''Main parent-offspring phasing processing chain within a family.'''
####################################################################################
def family_phaser(next_filter=None, debug=False, single_family=False, print_times=False):
    chain = [
             Filter(name='* Outer Duos', handle=__handle_outer_duos, debug=debug),
             # Phaser(name='* All Duos', handle=__handle_all_duos, debug=debug),
             FamilyPhaser(family_handle=__handle_parent_child, genotyped_families=True, min_children=0,
                          name='* Phased Parent->Child'),
             phased_samples_phaser  # In selected samples mode, reset all non-selected haplotypes back to the original values 
             ] 
    return new_phaser_chain(chain, name='Family Parent-Child', next_filter=next_filter,
                            debug=debug, print_times=print_times)

def family_child_comparison_phaser(next_filter=None, debug=False, print_times=False): 
    chain = [
             Filter(name='* Child Comparison', handle=__handle_child_comparison),
             phased_samples_phaser  # In selected samples mode, reset all non-selected haplotypes back to the original values 
             ]
    return new_phaser_chain(chain, name='Family Children', next_filter=next_filter,
                            debug=debug, print_times=print_times)
