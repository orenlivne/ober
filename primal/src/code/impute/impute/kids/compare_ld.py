#!/usr/bin/env python
'''
============================================================
Compare LD results from SNP and PLATO as a sanity check.

Created on July 15, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, numpy as np, matplotlib.pylab as P

path = os.environ['OBER_OUT'] + '/kids'
prefix = path + '/cytosnp/ld'
chrom = 22

# Load SNP, PLATO results, find common pairs and insert corresponding R^2 values into the x (SNAP)
# and y (PLATO) arrays.
snap = np.loadtxt(prefix + '/snap_results.txt', dtype=[('i', 'S12'), ('j', 'S12'), ('r', 'f4')], usecols=(0, 1, 3), skiprows=1)
plato = np.loadtxt(prefix + '/out.ld_calc_ldcalconly.txt', dtype=[('i', 'S12'), ('j', 'S12'), ('r', 'f4')], usecols=(1, 3, 4), skiprows=1)
s = dict(((x[0], x[1]), x[2]) for x in snap)
p = dict(((x[0], x[1]), x[2]) for x in plato)
x, y = zip(*[(s[k], p[k]) for k in set(s.keys()) & set(p.keys())])

# Plot
P.scatter(x, y)
P.xlabel('SNAP LD (CEU)')
P.ylabel('Plato LD (60 Hutt)')
P.title('Sanity Check of CytoSNP chr%d LD Using Two Sources' % (chrom,))
P.savefig(os.environ['OBER'] + '/doc/kids/ld_comparison.png')

