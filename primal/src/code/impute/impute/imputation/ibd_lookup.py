'''
============================================================
IBD dictionary lookup queries for SNP imputation.

Created on November 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
#import impute as im, numpy as np
from impute.ibd.segment import SegmentSet

#---------------------------------------------
# Constants
#---------------------------------------------        

#---------------------------------------------
# Methods
#---------------------------------------------
def segments_ibd_at_bp(ibd, bp, samples=None):
    '''Return the list of all samples that are reported as IBD in the closed base-pair segment [bp[0],bp[1]]. If a sample
    set ''samples'' is specified, accounts only for segments where those samples appear.'''             
    segments = ibd.segments_intersecting_bp_segment(bp)
    return SegmentSet(x for x in segments if samples & set(z[0] for z in x.samples)) \
        if samples else segments

def samples_ibd_at_bp(ibd, bp, samples=None):
    '''Return the list of all samples that are reported as IBD in the closed base-pair segment [bp[0],bp[1]]. If a sample
    set ''samples'' is specified, accounts only for segments where those samples appear.'''             
    segments = ibd.segments_intersecting_bp_segment(bp)
    return set(y[0] for x in segments for y in x.samples if samples & set(z[0] for z in x.samples)) \
        if samples else set(y[0] for x in segments for y in x.samples)

def samples_ibd_at_snp(ibd, snp, samples=None):
    '''Return the list of all samples that are reported as IBD in the closed base-pair segment [bp[0],bp[1]]. If a sample
    set ''samples'' is specified, accounts only for segments where those samples appear.'''             
    segments = ibd.segments_intersecting_snp(snp)
    return set(y[0] for x in segments for y in x.samples if samples & set(z[0] for z in x.samples)) \
        if samples else set(y[0] for x in segments for y in x.samples)

#---------------------------------------------
# Private Methods
#---------------------------------------------
