#!/usr/bin/env python
'''
============================================================
A data set of genotypes of a list of individuals in a
pedigree. This is a 3-D array (individual x SNP x allele).

Created on Jul 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
import matplotlib.pyplot as plt

x = np.arange(13, 0, -1)
y = 100*np.divide([741, 741, 741, 757, 785, 842, 875, 939, 982, 1032, 1072, 1123, 1150], 1416.)

plt.plot(x,y, 'bo-')
plt.xlabel('Min # Genotyped Children in Family')
plt.ylabel('% Nuclear Family members')

plt.savefig('family_stats.png')