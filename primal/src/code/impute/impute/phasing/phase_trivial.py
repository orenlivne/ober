'''
============================================================
Phase trivial cases in trios where at least one of the three
individuals (the proband or a parent) is homozygous at the
SNP of interest.

In some cases, missing data in a parent or child can be
imputed.

For shortcut, only filter handle() methods are defined here
and fed as a constructor argument to class Filter whenever
possible, i.e., whenever no additional state is required
in the Filter sub-class.

Created on July 26, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
from impute.tools import genotype_tools as gt, pedigree_tools as pt
from impute.data.constants import MISSING, ALLELES, CHILD, PATERNAL, MATERNAL
from impute.phasing.phase_core import new_phaser_chain
from chain import Filter
from impute.phasing.pre_processing import phased_samples_phaser

####################################################################################
def __handle_hom_entries(self, request):
    '''Phase all homozygous SNPs: find all homozygous SNPs in the genotype set
    and set their corresponding haplotypes to the genotypes.''' 
    problem = request.problem
    g = problem.genotype.data
    if request.params.selected_mode:
        # Phase only selected samples
        g = g[:, request.params.selected_samples, :]
        hom = np.where(gt.is_homozygous(g)[:, :])
        problem.haplotype.data[hom[0], request.params.selected_samples[hom[1]], :] = g[hom]
    else:
        # Phase all samples
        hom = gt.is_homozygous(g)[:, :]
        problem.haplotype.data[hom] = g[hom]
    return False

####################################################################################
def __handle_hom_parent(self, request):
    '''Homozygous parent => determine corresponding child haplotype entry.'''
    # Order matters: flag trio errors before duo treatment
    problem, params = request.problem, request.params
    __flag_trio_error(problem, params)
    for parent_type in ALLELES: __phase_duo_child(self, problem, params, parent_type)

def __flag_trio_error(problem, params):
    '''Flag incompatible parent-child trios as genotype errors. Note: there can still be
    duo errors, treated in Case A below.'''
    # Find genotyped children whose parents are genotyped
    # If selected samples, process only selected genotyped children with genotyped parents. 
    g, trio = problem.genotype.data, pt.selected_trios(problem, params)
    if not trio.size: return False
    gf, gm, gc = g[:, trio[:, PATERNAL], :], g[:, trio[:, MATERNAL], :], g[:, trio[:, CHILD], :]
    j = np.where(__equal_hom_parents_het_child(gf, gm, gc) | 
                 __hom_parent_not_in_hom_child(gf, gc) | 
                 __hom_parent_not_in_hom_child(gm, gc))
    # Flag errors in all trio members
    for i in xrange(3): problem.genotype_error(j[0], trio[j[1], i] , 'Same homozygous parents, het child')

def __phase_duo_child(self, problem, params, parent_type):
    '''A helper method that phases a child based on a single parent of type parent_type.'''
    # Find genotyped children whose parent is genotyped
    duo = pt.selected_duos(problem, params, parent_type)
    if not duo.size: return False
    g, h = problem.data
    gc, gp = g[:, duo[:, 0], :], g[:, duo[:, 1], :]
    # Restrict view to snps that are homozygous in parent  
    hom = np.where(gt.is_homozygous(gp)[:, :])
    parent_allele, gc_hom = gp[hom[0], hom[1], 0], gc[hom[0], hom[1], :]

    #------------------------------------------------------------------------------------  
    # Case A: Parent = (a,a), child = (b,b) and a != b (incompatible) ==> error
    #------------------------------------------------------------------------------------  
    j = np.where(np.logical_and(gt.is_homozygous(gc_hom)[:], gc_hom[:, 0] != parent_allele))[0]
    # Flag errors in both children and parents
    for i in xrange(2): problem.genotype_error(hom[0][j], duo[hom[1][j], i], 'Homozygous parent allele not found in child')
    
    #------------------------------------------------------------------------------------  
    # Case B: Parent = (a,a), child = (a,x) or (x,a), x in {0,1,2} (compatible) ==> 
    # set child parent hap to a and other hap to the other child genotype (x) 
    #------------------------------------------------------------------------------------
    # Note: h[array,array,:] is not a reference into h like h[scalar,scalar,:]. Thus,
    # when setting h, we must use h[original coordinates here] = ... . This occurs
    # several times in the code of this file.
    
    # Determine child haplotype corresponding to the hom parent 
    snps = gt.index_of(gc_hom, parent_allele)
    parent_allele_at_snps = parent_allele[snps]
    snp_index, child_index = hom[0][snps], duo[hom[1][snps], 0]
    h[snp_index, child_index, parent_type] = parent_allele_at_snps
    # Determine child haplotype corresponding to the other parent
    gc_rel = g[snp_index, child_index, :]
    other = np.where(gc_rel != np.transpose(np.tile(parent_allele_at_snps, (2, 1))))
    h[snp_index[other[0]], child_index[other[0]], 1 - parent_type] = gc_rel[other]
    
    #------------------------------------------------------------------------------------  
    # Case C: Parent = (a,a), child = (0,x) or (x,0) (potentially compatible) ==>
    # impute child to a and set child hap to a
    #------------------------------------------------------------------------------------
    for allele, snps in dict(zip(ALLELES, gt.index_first_missing(gc_hom))).iteritems():
        parent_value = parent_allele[snps]
        snp_original = hom[0][snps]
        gc[snp_original, hom[1][snps], allele] = parent_value
#        if self.debug:
#            print 'Imputing child', (snp_original, hom[1][snps], allele, parent_value)
#        problem.info.imputed_genotype.append((snp_original, hom[1][snps], allele, parent_value))
        h[snp_original, duo[hom[1][snps], 0], parent_type] = parent_value
    return False

'''Return indices in which Parent0 = Parent1 = (a,a) and child = (a,b), a != b (incompatible).'''
__equal_hom_parents_het_child = lambda gf, gm, gc: \
    gt.is_homozygous(gf)[:, :] & gt.is_homozygous(gm)[:, :] & (gf[:, :, 0] == gm[:, :, 0]) \
    & gt.is_heterozygous(gc)[:, :]
    
'''Return indices in which Parent = (a,b) and child = (b,b), a != b (incompatible).'''
__hom_parent_not_in_hom_child = lambda gp, gc: \
    gt.is_homozygous(gp)[:, :] & gt.is_homozygous(gc)[:, :] & (gp[:, :, 0] != gc[:, :, 0]) 

####################################################################################
def __handle_impute_parent(self, request):
    '''Child with two determined haps (a,b) and parent has (a,MISSING) or (MISSING,a) ==>
    impute parent to (a,b).
    '''
    # if request.params.selected_mode: return False
    for parent_type in ALLELES: __impute_duo_parent(request.problem, request.params, parent_type)

def __impute_duo_parent(problem, params, parent_type):
    '''Child with two determined haps (a,b) and parent has (a,MISSING) or (MISSING,a) ==>
    impute parent to (a,b). A helper method to process a single parent_type'''
    duo = pt.selected_duos(problem, params, parent_type)
    if not duo.size: return False
    g, h = problem.data
    gp, hc = g[:, duo[:, 1], :], h[:, duo[:, 0], :]

    # Restrict view to SNPs with full child haplotype
    child = np.where((hc[:, :, PATERNAL] != MISSING) & (hc[:, :, MATERNAL] != MISSING))
    gp_rel, hc = gp[child[0], child[1], :], hc[child[0], child[1], :]

    # Impute parent's first missing allele to be the corresponding child hap    
    for allele, snps in dict(zip(ALLELES, gt.index_first_missing(gp_rel))).iteritems():
        hc_value = hc[snps, :, parent_type]
#        if self.debug:
#            print 'Imputing parent', (full[0][snps], duo[full[1][snps],1], allele)
#        problem.info.imputed_genotype.append((full[0][snps], duo[full[1][snps],1], allele, hc_value))
        g[child[0][snps], duo[child[1][snps], 1], allele] = hc_value
    return False

def __handle_single_allele(self, request):
    '''Child with one determined hap and full data ==> determine other child allele.  
    Also, if parent of other allele has missing data, impute parent using child's data.'''
    problem, params = request.problem, request.params
    g, h, trio = problem.genotype.data, problem.haplotype.data, problem.kids_trios(params.selected_samples) if params.selected_mode else problem.trios()
    if not trio.size: return False
    child = trio[:, CHILD]
    gc = g[:, child, :]
    child_has_full_genotype = (gc[:, :, PATERNAL] != MISSING) & (gc[:, :, MATERNAL] != MISSING)
    
    for parent_type in ALLELES:
        other = 1 - parent_type
        # Other parent's genotype
        parent_other = trio[:, other]
        go = g[:, parent_other, :]
        # Corresponding child haplotypes
        hp, ho = h[:, child, parent_type], h[:, child, other]
        
        # Restrict view to relevant snps and samples:
        # child parent_type allele is determined; other allele is missing; child has full genotype  
        j = np.where(child_has_full_genotype & (hp != MISSING) & (ho == MISSING))
        gc_rel = gc[j[0], j[1], :]
        
        # Determine child haplotype corresponding to other parent's allele
        other_child_genotype = gc_rel[np.where(gc_rel != np.transpose(np.tile(hp[j[0], j[1], :], (2, 1))))]
        h[j[0], child[j[1]], other] = other_child_genotype

        # If parent of other allele has missing data, impute parent using child's data
        for allele, snps in dict(zip(ALLELES, gt.index_first_missing(go[j[0], j[1], :]))).iteritems():
            other_value = other_child_genotype[snps]
#            if self.debug:
#                print 'Imputing parent', (j[0][snps], parent_other[j[1]][snps], allele)
#            problem.info.imputed_genotype.append((j[0][snps], parent_other[j[1]][snps], allele,other_value))
            g[j[0][snps], parent_other[j[1]][snps], allele] = other_value
    return False

####################################################################################
'''Main trivial phasing processing chains - a facade.'''
def trivial_phaser(next_filter=None, debug=False, print_times=False):
    return new_phaser_chain([
                             Filter(name='* Hom SNPs', handle=__handle_hom_entries),
                             Filter(name='* Hom Parent', handle=__handle_hom_parent),
                             Filter(name='* Impute Parent', handle=__handle_impute_parent),
                             phased_samples_phaser  # In selected samples mode, reset all non-selected haplotypes back to the original values 
                             ], name='Trivial Cases', debug=debug, print_times=print_times)
