#!/usr/bin/env python
'''
============================================================
Estimate and plot the recombination rate lambda=lambda(f)
where f is the inbreeding coefficient. Per discussion with
Mark Abney on IBD HMM for haplotypes.  

Created on January 23, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib.pyplot as P, numpy as np, db_gene, impute as im

#---------------------------------------------
# Methods
#---------------------------------------------
def lambda_vs_f(id_coef=im.PhaseParam().id_coef_file):
    '''Return lambda as a discrete function of the inbreeding coefficient f for all Hutt samples.
    Outputs an array with f, lam columns.''' 
    # Read lambda, calculate f from Deltas
    k = __idcoef_dao(id_coef)
    a = np.array([np.concatenate(([y[0]], y[1])) for y in (k.id_coefs(x, x) for x in k.index)])
    inbreed = lambda d: 2 * (d[1] + 0.5 * (d[3] + d[5] + d[7]) + 0.25 * d[8]) - 1
    b = np.array([(inbreed(x), x[0]) for x in a])
    # Remove outliers with f~0
    b = b[np.where(b[:, 0] >= 1e-5)[0], :]
    return b

def lambda_mean(lam_array):
    '''Bin lambda, f and return the mean f, mean lambda and stddev of lambda in each bin.'''
    f, lam = lam_array[:, 0], lam_array[:, 1]
    bins = np.arange(0, max(f) * 1.1, 0.01)
    i = np.digitize(f, bins)
    c = range(1, len(bins))
    F = map(lambda m: 0.5 * (bins[m] + bins[m - 1]), c)
    L = map(lambda m: np.mean(lam[np.where(i == m)[0]]), c)
    S = map(lambda m: np.std(lam[np.where(i == m)[0]]), c)
    return F, L, S

def plot_lambda_vs_f(lam_array):
    '''Plot the recombination rate lambda vs. inbreeding coefficient.'''
    f, lam = lam_array[:, 0], lam_array[:, 1]
    F, L, S = lambda_mean(lam_array)
    P.clf()
    P.hold(True)
    P.plot(f, lam, 'bo')
    P.errorbar(F, L, yerr=S, fmt='ro-', linewidth=2, elinewidth=2)
    # P.ylim([0, 1.5])
    P.xlabel('Inbreeding Coefficient $f$')
    P.ylabel('Recombination Rate $\lambda$')
    # P.title('Recombination Rate vs. Inbreeding in the Hutterities')
    # P.show()

def plot_lambda_std(problem, id_coef=im.PhaseParam().id_coef_file):
    '''Plot lambda std dev vs. mean lambda in all children of all families in the problem object
    ''problem'.'''
    k = __idcoef_dao(id_coef)
    l = dict((x, k.id_coefs(x, x)[0]) for x in k.index)
    child_lam = [map(l.get, problem.pedigree.sample_id[np.array(list(y.children))]) for y in problem.families()]
    c = np.array([(np.mean(x), np.std(x)) for x in child_lam])
    P.clf()
    P.hold(True)
    P.scatter(c[:, 0], c[:, 1])
    P.xlabel('Sibs Mean $\lambda$')
    P.ylabel('Sibs Stddev $\lambda$')
    # P.show()

#---------------------------------------------
# Private Methods
#---------------------------------------------
__idcoef_dao = lambda id_coef: db_gene.snp.file_dao.IdCoefDao(id_coef)
