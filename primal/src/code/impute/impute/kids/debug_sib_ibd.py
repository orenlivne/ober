#!/usr/bin/env python
'''
============================================================
Compare CytoSNP and OmniExpress IBD segments for one sib.
Seems that CytoSNP finds segments and OmniExpress does not
which does not make sense since OmniExpress is more dense.

Phasing problem?

Created on July 31, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os
from impute.kids.hutt_kids import get_sib_ids, read_chip_problem

####################################################################################
if __name__ == '__main__':

    path = os.environ['OBER_OUT'] + '/kids'
    chrom = 22
    debug = 1
    num_processes = 1  # 4
    params = im.PhaseParam()
    
    sibs = get_sib_ids(path)
    sib = 800

    # Original Affy data
    print '*' * 80
    print 'Affy Original Phasing Output'
    print '*' * 80
    segments = im.examples.hutt_ibd_segments('hutt.phased.npz', 800, 1, 235, 0, debug=False)
    print segments
    
    for chip in ['affy', 'cytosnp']:  # , 'omniexpress']:
        print '*' * 80
        print 'Chip', chip
        print '*' * 80
        p = read_chip_problem(path, chrom, chip)    
        segments = im.problem_ibd_segments(p, 800, 1, 235, 0, debug=False)
        print segments
