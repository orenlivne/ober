#!/usr/bin/env python
'''
============================================================
Plot CGI imputation result - call rates. 

Created on Feb 17, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib.pyplot as P, os, impute as im, numpy as np

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    save_dir = os.environ['OBER'] + '/doc/imputation/cgi/v2'
    code_dir = os.environ['OBER_CODE'] + '/impute/batch/cgi'

    im.imputation.impute_stats.plot_cgi_call_rate(save_dir + '/count-all.txt',
                                                  save_dir + '/count-98.txt',
                                                  1415,
                                                  title='CGI Imputation Call Rates: All Hutterites Samples')
    P.savefig(save_dir + '/call-rate-all.png')

    # Selected samples with high call rates only
    im.imputation.impute_stats.plot_cgi_call_rate(save_dir + '/count-selected.txt',
                                                  save_dir + '/count-98.txt',
                                                  len(np.loadtxt(code_dir + '/selected.id', dtype=np.uint)),
                                                  title='CGI Imputation Call Rates: Selected Hutterites Samples')
    P.savefig(save_dir + '/call-rate-selected.png')
