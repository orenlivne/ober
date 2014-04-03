#!/usr/bin/env python
'''
============================================================
Hutterites vs. EUR AAF (Altternative Allele Frequenc)
plots for imputation paper. Uses CGI variant annotation data
file prepared as part of the annotation database.

Created on December 19, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pylab as P, os, numpy as np
from matplotlib.colors import LogNorm

# m = array of AAF values
# m[: ,0] = AAF EUR
# m[:, 1] = AAF CGI
# m[:, 2] = AAF imputed

def hist_novel_variants(m):
    '''Generate a plot of Hutt AAF for all variants with EUR AAF = 0.'''
    novel = (m[:, 0] < 0.00001)
    P.figure(1)
    P.clf()
    # P.hist(m[novel, 1], bins=20, log=True)
    P.hist(m[novel, 1:3], 10, normed=False, log=True, histtype='bar', label=['98 CGI', 'Imputed'])    
    P.xlabel(r'Hutterite AAF')
    P.ylabel('# Variants')
    P.title('Hutterite AAF for Non-EUR Variants (#variants = %d)' % (m.shape[0],))
    P.legend(loc='best')

def eur_vs_cgi_aaf(m):
    # Color proportional to #variants
    P.figure(2)
    P.clf()
    in_both = (m[:, 0] > 0) & (m[:, 1] > 0)
    P.hist2d(m[in_both, 1], m[in_both, 0], bins=40, norm=LogNorm())
    P.xlabel(r'Hutterite CGI AAF')
    P.ylabel('EUR AAF')
    P.title('EUR vs. Hutterite 98 CGI AAF (#variants = %d)' % (m.shape[0],))
    P.colorbar()

def eur_vs_imputed_aaf(m):
    # Color proportional to #variants
    P.figure(3)
    P.clf()
    in_both = (m[:, 0] > 0) & (m[:, 2] > 0)
    P.hist2d(m[in_both, 2], m[in_both, 0], bins=40, norm=LogNorm())
    P.xlabel(r'Hutterite Imputed AAF')
    P.ylabel('EUR AAF')
    P.title('EUR vs. Hutterite Imputed AAF (#variants = %d)' % (m.shape[0],))
    P.colorbar()

def eur_imputed_diff_aaf(m):
    '''A histogram of EUR AAF - Hutt imputed AAF.'''
    # Color proportional to #variants
    P.figure(4)
    P.clf()
    in_both = (m[:, 0] > 0) & (m[:, 2] > 0)
    P.hist(m[in_both, 0] - m[in_both, 2], log=True)
    P.xlabel('EUR AAF - Hutterites Imputed AAF ')
    P.ylabel(r'# Variants')

def plots_for_paper(m):
    doc_dir = os.environ['OBER'] + '/doc/paper_impute'
    hist_novel_variants(m)
    P.savefig(doc_dir + '/novel-hutt-aaf.png')
    eur_vs_cgi_aaf(m)
    P.savefig(doc_dir + '/eur-vs-cgi-aaf.png')
    eur_vs_imputed_aaf(m)
    P.savefig(doc_dir + '/eur-vs-imputed-aaf.png')
    eur_imputed_diff_aaf(m)
    P.savefig(doc_dir + '/eur-imputed-diff-aaf.png')

def load_data(): 
    m = np.loadtxt(os.environ['OBER_OUT'] + '/impute_cgi/annotations/aaf.txt', usecols=[1, 2, 3])
    # Sanity check: there should not be variants with CGI AAF ~ 0.5 and imputed AAF ~ 0
    m = m[(m[:,0] < 1) & (m[:,1] < 1)]
    novel = np.where(m[:, 0] < 0.00001)[0]
    j = ((m[novel, 1] > 0.45) & (m[novel, 2] == 0))
    print '# variants with large EUR AAF, small imputed AAF', len(np.where(j)[0]) 
    return m
    
####################################################################################
if __name__ == '__main__':
    pass
    m = load_data()
    plots_for_paper(m)
