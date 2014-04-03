#!/usr/bin/env python
'''
============================================================
Plot haplotype differences between sibs in the family 945
restricted problem.

Created on October 10, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib.pyplot as P
import impute as im

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    p = im.io.read_npz(im.itu.FAMILY945_ONE_PARENT_STAGE2).remove_nodes([2])

    f = p.first_family
    im.plots.plot_family_comparison(p, f, 1, template=1, children=[1, 2, 4])
    P.show()
    
#    P.figure(1)
#    im.plots.plot_all_diffs(p.haplotype, 1, 2, hap1_type=1, hap2_type=1)
#    P.show()
#    
#    P.figure(2)
#    im.plots.plot_all_diffs(p.haplotype, 1, 4, hap1_type=1, hap2_type=1)
#    P.show()
#    
#    P.figure(3)
#    im.plots.plot_all_diffs(p.haplotype, 2, 4, hap1_type=1, hap2_type=1)
#    P.show()
