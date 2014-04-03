#!/usr/bin/env python
'''
============================================================
Plot the recombination rate lambda=lambda(f) where f is
the inbreeding coefficient. Per discussion with Mark Abney
on IBD HMM for haplotypes.  

Created on January 23, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib.pyplot as P, impute as im, os

lam = im.hap_lambda.lambda_vs_f()
F, L, S = im.hap_lambda.lambda_mean(lam)

out_dir = os.environ['OBER'] + '/doc/ibd'

# Bin lambda, f and calculate mean and stddev of each bin so that we can see trends
P.close(1)
P.figure(1)
im.hap_lambda.plot_lambda_vs_f(lam)
# P.title('Recombination Rate vs. Inbreeding in the Hutterities')
P.show()
P.savefig(out_dir + '/lambda_vs_f.eps')

# Plot lambda std dev vs. mean lambda in family children
P.close(2)
problem = im.hutt('hutt.npz')
P.figure(2)
im.hap_lambda.plot_lambda_std(problem)
P.show()
P.savefig(out_dir + '/lambda_child.eps')
