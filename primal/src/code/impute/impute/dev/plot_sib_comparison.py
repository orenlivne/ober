#!/usr/bin/env python
'''
============================================================
Plot "correlation" matrices of haps of siblings in a family
with non-genotyped parents.   
 
Created on September 7, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, matplotlib.pyplot as plt, util
from impute import impute_test_util as itu
from impute.plot import plots
from impute.data import io

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    
    # Load nuclear family data
    p = io.read_npz(itu.FAMILY2003_STAGE3)
    f = p.families(genotyped=False)[0]
    old_printoptions = np.get_printoptions()
    np.set_printoptions(precision=2)
    
    children = np.array([x for x in f.children_list if p.is_genotyped(x)])
    h = p.haplotype.data
    plots.plot_hap_corr_matrix(h, children, (0, 0))
    plt.savefig('hap00.png')
    plots.plot_hap_corr_matrix(h, children, (0, 1))
    plt.savefig('hap01.png')
    
    util.set_printoptions(old_printoptions)
