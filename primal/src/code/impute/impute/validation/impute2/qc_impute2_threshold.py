#!/usr/bin/env python
'''
============================================================
Determine the optimal r^2 and MAF threshold - IMPUTE2 QC. 
This is a constrained optimization:

max call rate
s.t. accuracy <= eps (= 1% typically)

Usage: qc_impute2_threshold <file-name>
Use file-name = '-' to read from stdin.

Created on February 13, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import csv, sys

NUM_SAMPLES = 1415

# NUM_HETS_THRESHOLD = (0, 10, 20, 40, 80) 
MAF_THRESHOLD = (0, 30, 100) 

#MAF_THRESHOLD = (0, 0.005, 0.01, 0.02, 0.04)
#MAF_THRESHOLD = (0, 0.0107, 0.0367)

R2_THRESHOLD = (0.7, 0.8, 0.9)
'''Minimum concordance allowed'''
min_concordance = 0.99

'''Allocate an m x n matrix of zeros.'''
zeros = lambda m, n: [[0.0 for _ in xrange(n)] for _ in xrange(m)]

'''Convert het proportion to MAF (assuming HWE).'''
maf = lambda h: 0.5 * (1 - (1 - (2.0 * h) / NUM_SAMPLES) ** 0.5)

####################################################################################
if __name__ == '__main__':    
    try:        
        # Loop over file, accumulate call rate, concordance statistics
        m, n = len(MAF_THRESHOLD), len(R2_THRESHOLD)
        count = zeros(m, n)
        concordance = zeros(m, n)
        call_rate_increase = zeros(m, n)
        num_records = 0
        # Column numbers:
        # 9=#hets called by both imputation methods
        # 10=het_concordance
        # 11=call_rate_impute2
        # 13=impute2_info
        # 14=maf
        for x in csv.reader(sys.stdin if sys.argv[1] == '-' else open(sys.argv[1], 'rb'), delimiter='\t'):
            num_records += 1
            for i, maf_threshold in enumerate(MAF_THRESHOLD):
                if float(x[9]) >= maf_threshold:
                    r2 = float(x[13])
                    for j, r2_threshold in enumerate(R2_THRESHOLD):
                        if r2 >= r2_threshold:
                            count[i][j] += 1
                            concordance[i][j] += float(x[10])
                            call_rate_increase[i][j] += float(x[11])
                            # print (i, j), maf_threshold, x[9], x[10], x[11], x[13], x[14]
             
        # Report call rates vs. concordance results
        max_call_rate = 0
        i_max, j_max = -1, -1
        print '%-11s  %-7s  %-4s  %-4s' % ('Category', '# Vars', 'CR inc', 'Concord')
        print '-'*40
        for i, maf_threshold in enumerate(MAF_THRESHOLD):
            print '#hets >= %3d <==> MAF >=~ %.2f%%' % (maf_threshold, 100.*maf(maf_threshold),)
            for j, r2_threshold in enumerate(R2_THRESHOLD):
                cr = call_rate_increase[i][j] / num_records
                c = concordance[i][j] / count[i][j]
                print 'r2 >= %.2f   %-7d  %-.4f  %-.4f' % (r2_threshold, count[i][j], cr, c)
                if c >= min_concordance and cr > max_call_rate: i_max, j_max, max_call_rate = i, j, cr
            print ''
        
        # Print optimal solution
        i, j = i_max, j_max
        if i < 0: print 'No feasible solutions'
        else:
            cr = call_rate_increase[i][j] / num_records
            c = concordance[i][j] / count[i][j]
            print 'Max call rate with concordance >= %.2f: #hets >= %3d (MAF >~= %.2f%%), r2 >= %.2f' % (min_concordance, MAF_THRESHOLD[i], 100.*maf(MAF_THRESHOLD[i]), R2_THRESHOLD[j])
            print '%-11s  %-7d  %-.4f  %-.4f' % ('', count[i][j], cr, c)
            
    except (IOError, OSError):
        sys.exit(141)
