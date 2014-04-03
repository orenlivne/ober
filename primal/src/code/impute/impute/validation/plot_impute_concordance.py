#!/usr/bin/env python
'''
============================================================
Plot imputation SNP concordances for a validation data set.  

Created on June 5, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as P, os, numpy as np, sys, util

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    in_file = sys.argv[1]
    in_file_common = sys.argv[2]
    title = sys.argv[3]
    save_dir = sys.argv[4]  # os.environ['OBER'] + '/doc/imputation/rsng'
    fac = 0.99

    util.mkdir_if_not_exists(save_dir)
    data = np.loadtxt(in_file, usecols=[1, 2, 3, 4])
    concord = 1.0 - (1.0 * data[:, 3]) / (data[:, 1] - data[:, 2]) 

    if in_file_common != '-':
        data_common = np.loadtxt(in_file_common, usecols=[1, 2, 3, 4])
        concord_common = 1.0 - (1.0 * data_common[:, 3]) / (data_common[:, 1] - data_common[:, 2]) 

    snp_name = np.loadtxt(in_file, usecols=[0], dtype=str)
    
    P.figure(1)
    P.clf()
    P.hold(True)
    ax = P.subplot(111)
    ax.scatter(data[:, 0], concord, color='b', lw=0, label='All Samples')
    if in_file_common != '-':
        ax.scatter(data_common[:, 0], concord_common, color='r', lw=0, label='CGI Samples')
    ax.set_xscale('log')
    ax.set_xlim(1e-3, 1e0)
    P.draw()    
    P.xlabel('Minor Allele Frequency')
    P.ylabel('Concordance Rate')
    #P.ylim([min(concord) * fac, 1.01])
    P.title('Concordance Rate of %s Genotypes vs. Imputed\n#SNPs = %d, Mean %.2f%% Median %.2f%%' \
            % (title, data.shape[0], 100 * np.mean(concord), 100 * np.median(concord)))
    P.legend(loc='lower left', prop={'size': 10})
    # P.show()
    P.savefig(save_dir + '/%s-concordance.png' % (os.path.basename(in_file)))

    # Filter SNPs to high CGI-RS&G concordance
    if in_file_common != '-':
        good = np.where(concord_common >= 0.99)[0]
        P.figure(2)
        P.clf()
        P.hold(True)
        ax = P.subplot(111)
        ax.scatter(data[good, 0], concord[good], color='b', lw=0, label='All Samples')
        ax.scatter(data_common[good, 0], concord_common[good], color='r', lw=0, label='CGI Samples')
        ax.set_xscale('log')
        ax.set_xlim(1e-3, 1e0)
        P.draw()    
        P.xlabel('Minor Allele Frequency')
        P.ylabel('Concordance Rate')
        #P.ylim([min(concord_common[good] * fac), 1.001])
        P.title('Concordance Rate of %s Genotypes vs. Imputed\n#CGI-Concordant SNPs = %d, Mean %.2f%% Median %.2f%%' \
                % (title, len(good), 100 * np.mean(concord[good]), 100 * np.median(concord[good])))
        P.legend(loc='lower left', prop={'size': 10})
        # P.show()
        P.savefig(save_dir + '/%s-concordance-filtered.png' % (os.path.basename(in_file)))
        
        # Display non-(RS&G-CGI-concordant) SNPs
        data = np.concatenate((data, concord[np.newaxis].transpose()), axis=1)
        # i = np.argsort(data[:, 4])
        i = np.where(concord_common < 0.99)[0]
        stats = np.concatenate((i[np.newaxis].transpose(), data[i], concord_common[i][np.newaxis].transpose()), axis=1)  # [0:10]
        print '%-4s %-12s %-6s %-10s %-10s %-10s %-11s' % ('#', 'SNP', 'MAF', '#genotypes', '#discord', 'Imp.Concord' , 'CGI Concord')
        for row in stats:
            # print '%-4d %-12s %-6.3f %-10d %-10d %-10.3f' % (row[0], snp_name[i[row[0]]], row[1], row[2] - row[3], row[4], row[5])
            print '%-4d %-12s %-6.3f %-10d %-10d %-10.3f %-11.3f' % (row[0], snp_name[row[0]], row[1], row[2] - row[3], row[4], row[5], row[6])
