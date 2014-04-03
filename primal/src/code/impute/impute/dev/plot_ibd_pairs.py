#!/usr/bin/env python
'''
============================================================
Plot statistics on the number of Hutt pairs that are IBD
over entire genetic regions. This can guide IBD segment set
pruning and optimization.   
 
Created on March 11, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, matplotlib.pyplot as P, impute as im, os

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    out_dir = os.environ['OBER'] + '/doc/imputation/validation/ibd-optimization'
    
    # Data format: [region start snp, region stop snp, # IBD pairs in 1415x1415 IBD set]
    # Data obtained manually with commands like
    #i=50; j=3168; cat segments.out | awk -v i=$i -vj=$j '($1 <= i) && (j <= $2)' | wc -l
    # (or into file and then wc -l)

    pairs = [((50, 3168), 1428), ((550, 2650), 7554), ((1050, 2150), 34630), ((1300, 1900), 62017), ((1425, 1775), 94152)]
    
    p = im.hutt('hutt.npz')
    cm = p.info.snp['dist_cm']
    chrom = p.info.snp['chrom'][0]
    
    l = [cm[x[0][1]] - cm[x[0][0]] for x in pairs]
    num_pairs = [x[1] for x in pairs]    
    
    P.figure(1)
    P.clf()
    P.semilogy(l, num_pairs, 'bo-')
    P.grid(True)
    P.title('Chromosome %d' % (chrom,))
    P.xlabel('Region Length [cM]')
    P.ylabel('# IBD Pairs')
    P.show()
    P.savefig('%s/num_ibd_pairs.chr%d.png' % (out_dir, chrom,))
    
    # Analyze a specific segment set for 15cM region in the middle of chr22 (34630 segments)
    # awk {'printf "%d %d %d %d %d %d %d %d\n", 0, 3217, 16484792, 51156933, $5,$6,$7,$8'} segments.15cm.out  > segments.15cm.fake.out
    a = im.segment.SegmentSet.load(open(os.environ['OBER_OUT'] + '/ibd/chr22/segments.15cm.fake.out', 'rb'))
    a.group_to_disjoint()
    G = np.array([len(x.samples) for x in a])
    P.figure(2)
    P.clf()
    P.hist(G, max(G))
    P.title('Haplotype Groups, Chromosome %d (#groups=%d, #haps=%d)' % (chrom, len(G), sum(G)))
    P.xlabel('Haplotype Group Size')
    P.ylabel('Frequency')
    P.show()
    P.savefig('%s/hap_groups.chr%d.15cm.png' % (out_dir, chrom,))
    
