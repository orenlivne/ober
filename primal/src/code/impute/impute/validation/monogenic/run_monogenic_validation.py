#!/usr/bin/env python
'''
============================================================
Run manual imputation on monogenic SNPs (Jessica's rare
SNPs) for which there were discordances between the CGI and
Daokta genotypes. 

Created on June 12, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, numpy as np, impute as im

im.v.monogenic.monogenic_validation.pipeline_monogenic_validation(work_dir=os.environ['OBER_OUT'] + '/requests/monogenic/work',
                                                                  index_segments_dir=os.environ['OBER_OUT'] + '/requests/monogenic/work/index_segments',
                                                                  regenerate_segments=False,  # True,
                                                                  #snps=np.array([6]),
                                                                  debug=1,
                                                                  debug_sample=944)
