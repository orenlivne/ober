#!/usr/bin/env python
'''
============================================================
Compare the GERMLINE IBD segments on top of the nuclear
family recombination picture we calculate with the
ibd_family module algorithm.
 
Created on October 26, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im
#import matplotlib.pyplot as plt

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    p = im.io.read_npz(im.itu.FAMILY4_STAGE3)
    m = im.ig.ibd_germline(p, p.first_family.member_list)
    m.pprint_segments(True)
    d = m.group_by_segment()
    print d
    
    # Number of distinct haplotypes found
    num_haps = max(len(v) for v in d.itervalues())
    
    # Compute and plot differences    
#    plt.figure(1)
#    plots.plot_all_diffs(h, a, b, hap2_type=constants.MATERNAL)
#    plt.savefig('hap_diff_%d_%d.png' % (a, b, ))
