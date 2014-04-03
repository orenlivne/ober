#!/usr/bin/env python
'''
============================================================
Variant count classes, business logic and plots.

Based on http://matplotlib.org/examples/pylab_examples/bar_stacked.html

Created on April 2, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, os
import matplotlib.pyplot as plt
from annotate.hutt_dao import VariantDao, VariantSummaryReport

#---------------------------------------------
# Constants
#---------------------------------------------
# colors = ['b', 'g', 'r', 'c', 'm', 'y', 'g', 'b', 'k', 'r', 'c']  # Colors of stacked groups
COLOR_MAP = plt.get_cmap('gist_rainbow')

#---------------------------------------------
# Methods
#---------------------------------------------
def variant_summary_bar_chart(report):
    '''Generate a bar chart from a variant summary report object (VariantSummaryReport).'''
    labels, maf, count = report._categories, report._maf, report._count
    num_groups, N = count.shape
    width = 1  # the width of the bars: can also be len(x) sequence
    xticks = width * (np.arange(len(maf) - 1) + 0.5)
    
    fig = plt.figure(num=1, figsize=(9.5, 4.1), dpi=150)  # , facecolor='w', edgecolor='k')
    fig.clf()
    ax = plt.subplot(111)
    # colors = [COLOR_MAP(1.*i / num_groups) for i in xrange(num_groups)]
    colors = ['r', 'y', 'g', 'b', 'm', 'gray']
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    p = [ax.bar(xticks, count[:, group],
                bottom=(0 if group == 0 else np.sum(count[:, :group], axis=1)),
                log=True, width=0.5 * width, color=colors[group]) 
         for group in xrange(N)]
    
    plt.xlabel('Minor Allele Frequency', fontsize=14)
    plt.ylabel('# Variants', fontsize=14)
    plt.xticks(xticks + 0.25 * width, [('<%.2f' % (maf[i + 1]) if i == 0 else '%.2f-%.2f' % (maf[i], maf[i + 1])) for i in xrange(len(maf) - 1)], fontsize=7)
    # plt.yticks(np.linspace(1, np.max(count), 10), fontsize=18)
    
    # Put a legend to the right of the current axis
    ax.legend(tuple(x[0] for x in p), labels, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    plt.xlim([-width * 0.05, max(xticks) + width * 1.05])
    
def generate_reports_for_paper(db_url='mysql://hutt:hutt@127.0.0.1/hutt', save_dir=os.environ['OBER'] + '/doc/paper_impute/supplementary'):
    '''Generate all variant reports (data file + bar charts) for imputation paper using
    the database whose connection is specified by the URL ''db_url''.''' 

    # Load data from file
    for group in ('region', 'coding'):
        variant_dao = VariantDao(db_url)
#        report = variant_dao.variant_count_report(group)
#        report.save(open('%s/variant-count-%s.dat' % (save_dir, group), 'wb'))
        report = VariantSummaryReport.load(open('%s/variant-count-%s.dat' % (save_dir, group), 'rb'))
        # Generate figure
        variant_summary_bar_chart(report)
        # Save figure to file
        for extension in ('png', 'eps'):
            plt.savefig('%s/variant-count-%s.%s' % (save_dir, group, extension))
