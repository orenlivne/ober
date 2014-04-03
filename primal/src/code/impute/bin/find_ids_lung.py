#!/usr/bin/env python
'''
============================================================
Find Lung function study ID sublist of the entire Hutterites
problem set ID list. Output indices into the latter.

Created on Feb 18, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, numpy as np, impute as im
#from optparse import OptionParser

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    d = os.environ['OBER_DATA'] + '/lung'
    s = im.io_pedigree.read(d + '/hutt-lung.pdg.tfam', genotyped_id_file=d + '/hutt-lung-samples.txt')
    p = im.hutt('hutt.stage5.npz')
    i = np.array([i in s.sample_id[0:s.num_genotyped] for i in p.pedigree.sample_id[0:p.pedigree.num_genotyped]])
    np.savetxt(sys.stdout, np.where(i)[0], fmt='%d')
