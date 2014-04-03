#!/usr/bin/env python
'''
============================================================
Determining thresholds during QC. 
Fit a curve to the concordance.

Created on November 1, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, os
import matplotlib.pylab as P 
from scipy.optimize.minpack import curve_fit, fsolve

VARIANT_NAME = {'snp': 'SNVs', 'del': 'Deletions', 'ins': 'Insertions', 'sub': 'Substitutions'}
TYPE_NAME = {'sgl': 'Singleton', 'nc': ''}

'''Concordance vs. call rate model.'''
# model = lambda x, a, b, c: a * np.exp(-b * x) + c
# model = lambda x, a, b: a + b / (1 - x)  # + c / (1 - x) ** 2  # a * np.exp(-b * x) + c
# model = lambda x, a, b: a + b * x ** 0.5  # + c / (1 - x) ** 2  # a * np.exp(-b * x) + c
model = lambda x, a, b, c: a + b * x ** c  # + c / (1 - x) ** 2  # a * np.exp(-b * x) + c

def fit_curve(r, n, c):
    '''Fit the data c(r) where c=concordance and r=call rate using an exponential approximation f.
    The fit is weighted using a normal approximation to a sum of binomial distributions, each representing
    a single events of a SNP being concordant for a single pair of samples at a single IBD2 segment.'''
    # print r, c
    (p, _) = curve_fit(model, r, c)  # , sigma=np.maximum(1e-5, np.sqrt(c * (1 - c) / n)))
    return model(r, *p), p

def plot_data_and_model(r, c, f, cutoff, title):
    '''Generate a plot of the data c(r) and fitted curve f.'''
    P.figure()
    P.clf()
    P.hold(True)
    P.plot(r, 1 - c, '.', label='Data', markersize=10)
    # P.plot(r, f, 'r-', label='Fit', linewidth=2)
    # P.axvline(x=cutoff, linewidth=2, color='k')
    # P.xlabel(r'$\beta$')
    P.xlabel('Call Rate', fontsize=18)
    # P.ylabel('Frequency')
    P.ylabel('Error Rate', fontsize=18)
    # P.title('%s: cutoff at %.2f' % (title, cutoff))
    P.title('%s' % (title,), fontsize=18)
    # P.legend(loc=0)
    # P.plot(x[knee], y[knee], 'ko', markersize=8)

def fit_data(file_name, usecols, title, out_name, threshold=0.99, out_dir=os.environ['OBER'] + '/doc/qc'):
    '''Fit a curve QC data and calculate a cutoff.'''

    # Load data
    data = np.loadtxt(file_name, skiprows=1, usecols=usecols, dtype=float)
    r, n, c = data[::2, 0], data[::2, 1], data[::2, 2]
    has_data = (n > 0)
    r, n, c = r[has_data], n[has_data], c[has_data]
    r = 1.0 - r / r[-1]  # Convert no call count to call rate

    # Fit model
    f, p = fit_curve(r, n, c)
    g = lambda x: model(x, *p) - threshold
    cutoff = fsolve(g, 0.5)
    
    # Generate plot
    plot_data_and_model(r, c, f, cutoff, title)
    P.savefig('%s/%s.png' % (out_dir, out_name))
    
    if cutoff >= 1:
        print '%-40s: don\'t use at all' % (title,)
    else:
        print '%-40s: %.2f   mean %.3f' % (title, cutoff, 1 - np.mean(c))
    return r, n, c

####################################################################################
if __name__ == '__main__':
    in_dir = os.environ['OBER_OUT'] + '/qc/datafiles'
    threshold = 0.99
    print 'Call Rate Cutoffs and Mean Error Rates for concordance threshold = %.3f' % (threshold,)
    P.close('all')
    old_settings = np.seterr(all='ignore')
    for t in ('nc', 'sgl'):
        type_name = TYPE_NAME[t]
        print t
        for variant_type in ('snp', 'del', 'ins'):#, 'sub'):
            variant_name = VARIANT_NAME[variant_type]
            file_name = in_dir + '/HH-fltr_chr000_%s_nc%s.tsv' % (variant_type, '' if t == 'nc' else '_sgl')
            fit_data(file_name, [0, 3, 4], '%s' % (variant_name,), '%s_%s_rs' % (variant_type, t), threshold=threshold)
            fit_data(file_name, [0, 7, 8], '%s' % (variant_name,), '%s_%s_no_rs' % (variant_type, t), threshold=threshold)
#             fit_data(file_name, [0, 3, 4], '%s %s with RS Numbers' % (type_name, variant_name,), '%s_%s_rs' % (variant_type, t), threshold=threshold)
#             fit_data(file_name, [0, 7, 8], '%s %s without RS Numbers' % (type_name, variant_name,), '%s_%s_no_rs' % (variant_type, t), threshold=threshold)
        print ''
