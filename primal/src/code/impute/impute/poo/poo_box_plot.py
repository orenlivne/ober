#!/usr/bin/env python
'''
============================================================
Generate some box plots for POO-GWAS hits. 

Created on January 22, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import csv, numpy as np, os, sys, itertools as it, matplotlib.pyplot as P
# from optparse import OptionParser

#---------------------------------------------
# Constants
#---------------------------------------------
'''# of metadata columns in an itabix line.'''
NUM_METADATA_COLS = 8 

#---------------------------------------------
# Methods
#---------------------------------------------
def hit_boxplot(path, hit_name):
    metadata = list(csv.reader(open('%s/%s.metadata' % (path, hit_name,), 'rb')))[0]
    x = np.loadtxt('%s/po-%s.txt' % (path, hit_name,),
                   dtype=[('id', 'i4'), ('p_allele', 'S1'), ('m_allele', 'S1'), ('phenotype', 'f8')])
    g = list(it.product(np.setdiff1d(np.unique(x['p_allele']), ['0']),
                        np.setdiff1d(np.unique(x['m_allele']), ['0'])))
    categories = [x[np.where((x['p_allele'] == a1) & (x['m_allele'] == a2))]['phenotype'] for a1, a2 in g]
    genotypes = map(''.join, g)
    P.figure()
    P.clf()
    P.boxplot(categories)
    
    
    # Set labels to genotypes + number of samples in each genotype group
    labels = map(lambda (x, y): '%s\n(%d)' % (x, y), zip(genotypes, map(len, categories)))
    P.xticks(np.arange(1, len(genotypes) + 1), labels)
    P.xlabel('Genotype')
    P.ylabel(metadata[1])
    P.title('%s (chr%d:%d)\n%s, %s' % (metadata[13], int(metadata[2]), int(metadata[3]), metadata[17], metadata[18]))
    P.savefig('%s/%s.png' % (path, hit_name))
    
####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    path = sys.argv[1]
    for hit_name in sys.argv[2].split():  # ('eosinophil', 'LDL'):
        hit_boxplot(path, hit_name)
