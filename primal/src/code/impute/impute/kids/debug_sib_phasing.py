#!/usr/bin/env python
'''
============================================================
Debug why CytoSNP re-phasing of a sample phased it wrong
already at stage 1.

The segments obtained on different platforms (chips) should
be similar (in their base-pair beginning and end positions). 

Created on August 2, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os, numpy as np, matplotlib.pyplot as P, sys
from impute.kids.hutt_kids import get_sib_ids, read_chip_problem, chip_data_set
from impute.phasing.phase_trivial import trivial_phaser

####################################################################################
if __name__ == '__main__':

    path = os.environ['OBER_OUT'] + '/kids'
    chrom = 22
    debug = 1
    num_processes = 1  # 4
    params = im.PhaseParam()
    
    sibs = get_sib_ids(path)
    sib = 842  # Has high imputation call rate but high re-phasing error rate

    chip = 'cytosnp'
    p = im.io.read_npz('%s/%s/chr%d/%s.npz' % (path, chip, chrom, chip_data_set(chip)))
    
    # Run stage 1 so that we can step through the code and find why SNP 1700
    # is phased to 2,1 instead of 1,2
    phaser = trivial_phaser()
    problem = phaser.run(p, im.PhaseParam(selected_samples=np.array([842]), debug=True, print_times=True))
