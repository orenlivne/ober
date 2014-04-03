#!/usr/bin/env python
'''
============================================================
Plot haplotype differences between two individuals in the
family test data set.
 
Created on Jul 18, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
#import numpy as np
import matplotlib.pyplot as plt
from impute import impute_test_util as itu
from impute.plot import plots
from impute.data import constants

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    
    # Input arguments
    h = itu.Templates.problem_family(itu.FAMILY7).haplotype
    a = 3
    b = 1
    
    # Compute and plot differences    
    plt.figure(1)
    plots.plot_all_diffs(h, a, b, hap2_type=constants.MATERNAL)
    plt.savefig('hap_diff_%d_%d.png' % (a, b, ))
