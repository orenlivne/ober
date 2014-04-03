#!/usr/bin/env python
'''
============================================================
A driver that generates all variant summary reports (by 
functional annotation) for the imputation paper.

Created on April 2, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib
matplotlib.use('Agg')

from annotate.variant_count import generate_reports_for_paper

####################################################################################
if __name__ == '__main__':
    generate_reports_for_paper(db_url='mysql://hutt:hutt@127.0.0.1/hutt')
