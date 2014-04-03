#!/usr/bin/env python
'''
============================================================
Debug POO on an affy SNP that was reported to have a large
transmission distortion among non-QF samples.

Created on January 27, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os  # , sys
from collections import Counter

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    snps = [409, 746]  # Affy snps exhibiting transmission distortion
    p, t = im.v.iv.impute_chrom(22, snps=snps, debug=1, input_dir=os.environ['OBER_OUT'] + '/phasing/chr22')
    # im.cgi.io_cgi.write_imputed(t, sys.stdout, poo_phase=p.haplotype.poo_phase)
    F = p.families(genotyped=True)
    GENOTYPES = ((1, 1), (1, 2), (2, 2))
    for snp in xrange(len(snps)):
        print 'SNP #%d' % (snp,)
        print '%-7s %-7s %-7s %-7s %-7s %-7s' % ('G(F)', 'G(M)', '#Fam', '#Kids', '12-Kids', '21-Kids')
        for gf in GENOTYPES:
            for gm in GENOTYPES:
                c = [(Counter(tuple(x) for x in p.h[snp, f.children_array]), f) for f in F 
                     if all(p.g[snp, f.father] == gf) and all(p.g[snp, f.mother] == gm)]
                num_children = sum(sum(v for k, v in x[0].iteritems() if 0 not in k) for x in c)
                print '%-7s %-7s %-7d %-7d %7.2f %7.2f' % \
                (repr(gf), repr(gm), len(c), num_children, \
                 sum(x[0][(1, 2)] for x in c) / float(num_children), \
                 sum(x[0][(2, 1)] for x in c) / float(num_children))
        c = Counter(tuple(x) for f in F for x in p.h[snp, f.children_array])
        print '%-15s %-7d %-7d %-7d %-7d' % ('TOTAL', len(F), sum(c.values()), c[(1, 2)], c[(2, 1)])
