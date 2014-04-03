#!/usr/bin/env python
'''
============================================================
Validation on non-QC affy SNPs.

Created on Oct 19, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
# For headless environments, must be imported before everything else
import matplotlib
matplotlib.use("Agg")
import sys, impute as im, os

####################################################################################
if __name__ == '__main__':
    im.v.iv.pipeline_validation_hutt_plink(sys.argv[1], os.path.dirname(sys.argv[1]), save_location=sys.argv[2])
