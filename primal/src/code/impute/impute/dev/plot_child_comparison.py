#!/usr/bin/env python
'''
============================================================
Example of recombination location identification and
phaing within a nuclear family in main scope.
 
Created on August 10, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
from impute import impute_test_util as itu
from impute.plot import plots
import impute.ibd.ibd_child as ic
from impute.data import constants, io

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    
    # Load nuclear family data
    p = io.read_npz(itu.FAMILY13_STAGE2)
    f = p.families()[0]
    
    # Calculate recombinations
    parent_type = constants.PATERNAL
    template = 2
    parent = f.parents[parent_type]
    (_, _, info) = ic.child_recombinations(p, f, parent_type, template, remove_errors=True)

    # Different ways to visualize the results
    print info.recombination_snp
    info.plot(template=True)
    info.plot(template=False)
    plots.plot_family_comparison(p, f, parent_type, template=template)
    plots.plot_family_comparison(p, f, parent_type)
