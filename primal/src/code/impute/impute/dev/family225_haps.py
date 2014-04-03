#!/usr/bin/env python
'''
============================================================
Display haplotypes generated in the nuclear family (225,81).
Debuggin stage 3 - Children comparison

Created on September 11, 2012 - my 10-year anniversary in the US!
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, util, numpy as np
from impute.data.constants import PATERNAL, ALLELES, MATERNAL
from impute.phasing.phase_family import family_child_comparison_phaser

def print_haps(h, hh, parent_type, snps=None):
    '''Print the parent haplotype pair and corresponding child (single) haplotypes for a SNP range.'''
    snps = snps if snps is not None else np.arange(0, p.num_snps) 
    parent = f.parents[parent_type]
    data = np.concatenate((snps[:, np.newaxis], h [snps, parent, :], h [snps, :, parent_type][:, children],
                           snps[:, np.newaxis], hh[snps, parent, :], hh[snps, :, parent_type][:, children]),
                          axis=1)
    print 'Haplotypes, parent type', parent_type
    print data

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    old_printoptions = np.get_printoptions()
    np.set_printoptions(threshold=np.nan, linewidth=120)
    
    # Load data
    # Before: p, h
    # After: q, hh
    p = im.io.read_npz(im.itu.FAMILY225_STAGE1)
    q = im.io.read_npz(im.itu.FAMILY225_STAGE1)
    phaser = family_child_comparison_phaser(debug=True)
    (g, h) = p.data
    (gg, hh) = q.data
    f = p.families(genotyped=True)[0]
    (father, mother) = f.parents
    children = np.array(f.children_list)
    
    phaser.run(q)
    
    # Print a portion of the father and corresponding children haplotypes
    snps = np.arange(2094, 2118) #np.arange(0, p.num_snps)
    for parent_type in ALLELES:
        print_haps(h, hh, parent_type, snps)
    
    # Recombinations
    comparator = im.ic.ChildComparator(p, f)
    template = {PATERNAL: 5, MATERNAL: 3}
    for parent_type in ALLELES:
        (_, _, info) = comparator.child_recombinations(parent_type, template[parent_type], remove_errors=False)
        print 'Recombinations, parent_type', parent_type
        print info.recombination_snp
    
    # Plot recombinations
    #im.plots.plot_all_family_comparisons(q, f) 

    util.set_printoptions(old_printoptions)
    
