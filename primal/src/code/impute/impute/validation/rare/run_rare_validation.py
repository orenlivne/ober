#!/usr/bin/env python
'''
============================================================
Validation on non-QC affy SNPs.

Created on Oct 19, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
# For headless environments
import matplotlib
matplotlib.use("Agg")

import sys, impute as im, os

####################################################################################
if __name__ == '__main__':
    # Input arguments
#    data = os.environ['OBER_OUT'] + '/requests/rare-combined/rare-combined'
    data = sys.argv[1]  # Name of data set 
    out = os.path.dirname(data)
    save_location = sys.argv[2]  # Output directory to save plots under

    im.v.iv.pipeline_validation_hutt_plink(data, out, snp_style='discrete', save_location=save_location, debug=True)
