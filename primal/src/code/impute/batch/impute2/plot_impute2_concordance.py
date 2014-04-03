#!/usr/bin/env python
'''
============================================================
Plot IMPUTE2 SNP concordance with our imputation.

Parameters: <impute2-stats-file-prefix> [dir-to-save-plots-in]

Created on September 17, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, impute.validation.impute2.impute2_validation as v

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    in_file = sys.argv[1]
    save_dir = None if len(sys.argv) < 3 else sys.argv[2]
    
    stats = v.impute2_concordance(in_file)
    # Print statistics
    # print 'Mean concordance: all genotypes %.2f het genotypes %.2f' % (stats[0], stats[1])
    v.plot_impute2_concordance(stats, save_dir=save_dir)
    
