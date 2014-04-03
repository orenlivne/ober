#!/usr/bin/env python
'''
============================================================
Intersect list of individuals with missing imputed genotypes
with the list of people that requested results back.
Requested by Rebecca Anderson. 

Created on March 13, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
path = 'c:/users/oren/Desktop'
# Format of imputed missing imputed genotypes: snp name, findiv
m = np.loadtxt(path + '/missing.txt', dtype=[('snp', 'S12'), ('id', 'i8')], delimiter=',')
# Assuming a header line followed by a list of findivs, one per line
r = set(np.loadtxt(path + '/rebecca.csv', dtype=int, skiprows=1))
# Intersect and save result back as CSV
np.savetxt(path + '/missing-and-reported-to.txt',
           np.array([x for x in m if x['id'] in r]), delimiter=',', fmt='%s')
