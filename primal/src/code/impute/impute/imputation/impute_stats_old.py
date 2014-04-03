#!/usr/bin/env python
'''
============================================================
Calculate imputation statistics; generate summary plots.

Created on December 3, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib.pyplot as P, numpy as np, util, impute as im, collections
from impute.tools import recode

#---------------------------------------------
# Methods - Rare Allele Imputation Stats
# (Older code)
#---------------------------------------------
def concordance(g1, g2, data_filter=im.imputation.istat.CALLED, samples=None):
    '''Calculate the concordance (# of filled entries that agree) between two genotype data arrays g1,g2.
    Returns an array of SNP indices that have data for comparison, and the corresponding concordance rates.'''
    if samples is not None: g1, g2 = g1[:, samples, :], g2[:, samples, :]
    r1, r2, groups = recode_single(g1, g2, data_filter)
    # print np.concatenate((g1[groups], g2[groups]), axis=1)
    groups = groups[0]
    concordant, snps = util.sum_by_group((r1 == r2).astype(np.byte), groups)
    discordant, _ = util.sum_by_group((r1 != r2).astype(np.byte), groups)
    c = np.zeros((max(groups) + 1,), dtype=np.float)
    c[snps] = 1.0 * concordant / (concordant + discordant)
#    return snps, 1.0 * concordant / (concordant + discordant)
    return c

def all_concordances(genotypes, imputed):
    '''Compute all measures of interest between the genotyped (iPlex) and imputed data that uses CGI as a
    training set.'''
    g, i = genotypes.data, imputed.imputed_data
    iPlex, cgi = (g[:, imputed.sample_index, :], imputed.genotype.data)
    maf = 100.0 * im.gt.allele_frequencies_by_snp(g)[1, :]

    c_training = concordance(iPlex, cgi, data_filter=im.imputation.istat.CALLED)
    c_training_rare = concordance(iPlex, cgi, data_filter=im.imputation.istat.HAS_MINOR_ALLELE)
    
    c_imputed = concordance(g, i, data_filter=im.imputation.istat.CALLED, samples=imputed.sample_index_to_impute)
    c_rare = concordance(g, i, data_filter=im.imputation.istat.CALLED_AND_HAS_MINOR_ALLELE, samples=imputed.sample_index_to_impute)
    
    return maf, c_training, c_training_rare, c_imputed, c_rare

def recode_single(g1, g2, data_filter):
    '''Recode a genotype pair to a single genotype code as in the recode module. Restrict to entries
    where both are called.'''
    r1, r2 = recode.recode_single_genotype(g1), recode.recode_single_genotype(g2)
    called = np.where(data_filter(r1, r2))
    return r1[called], r2[called], called

def num_training_hets(g, sample_index):
    '''Return an array with the # hets among training sample set ''sample_index'' for each of the SNPs
    in the g data array.'''
    x = collections.Counter(np.where(im.recode.recode_single_genotype(g.data[:, sample_index, :]) == 3)[0])
    y = np.zeros((max(x.iterkeys()) + 1,), dtype=np.uint32);
    for k, v in x.iteritems():
        y[k] = v
    return y

def snp_comparison_data(g, t, s):
    '''Debugging genotype vs. imputation data for SNP index s, stacked as a 2-D array.''' 
    a = np.concatenate((g.data[s, :, :], t.imputed_data[s, :, :], t.pedigree.genotyped_sample_id()[np.newaxis].transpose(), t.pedigree.genotyped_sample_index()[np.newaxis].transpose()), axis=1)
    # Hets
    hets = a[(a[:, 0] > 0) & (a[:, 1] > 0) & (a[:, 2] > 0) & (a[:, 3] > 0) & ((a[:, 0] + a[:, 1] == 3) | (a[:, 2] + a[:, 3] == 3)), :]
    R = t.sample_index
    training = np.concatenate((g.data[s, R, :], t.imputed_data[s, R, :], t.pedigree.genotyped_sample_id()[R][np.newaxis].transpose(), R[np.newaxis].transpose()), axis=1)
    return a, hets, training

#---------------------------------------------
# Methods - Plots
#---------------------------------------------
def plot_call_rate(stats):
    '''Plot imputation call rate vs. minor allele frequency).'''
    P.clf()
    P.scatter(np.min(stats['frequency'], axis=1), stats['call_rate'])
    P.title('Imputation Call Rate vs. MAF');
    P.xlabel('MAF')
    P.xlim([0, 0.5])
    P.ylim([P.ylim()[0], 1.0])
    P.ylabel('Call Rate')

def plot_all_concordances_vs_maf(snp_labels, all_concordances, threshold=0.0):
    '''Plot several concordance measured between measured genotypes and imputed genotypes.'''
    maf, c_training, c_training_rare, c_imputed, c_rare = all_concordances
    plot_concordance(snp_labels, maf, c_training, 'b', 'WGS Samples')
    plot_concordance(snp_labels, maf + 0.01, c_training_rare, 'k', 'WGS Samples, Rare')
    plot_concordance(snp_labels, maf + 0.02, c_imputed, 'g', 'Imputed Samples', threshold=threshold)
    plot_concordance(snp_labels, maf + 0.03, c_rare, 'r', 'Imputed Samples, Rare', threshold=threshold)

#    P.ylim([-0.05, 1.05])
    P.xlabel('Minor Allele Frequency [%]')
    P.ylabel('iPlex-Imputation Concordance')
    P.title('Imputation Accuracy for Rare Variants in the Hutterites')
    P.legend(loc='lower right', prop={'size': 10})

def plot_all_concordances_vs_num_hets(snp_labels, num_hets, all_concordances, threshold=0.0):
    '''Plot several concordance measured between measured genotypes and imputed genotypes.'''    
    _, c_training, c_training_rare, c_imputed, c_rare = all_concordances
    plot_concordance(snp_labels, num_hets, c_training, 'b', 'WGS Samples')
    plot_concordance(snp_labels, num_hets + 0.005, c_training_rare, 'k', 'WGS Samples, Rare')
    plot_concordance(snp_labels, num_hets + 0.01, c_imputed, 'g', 'Imputed Samples', threshold=threshold)
    plot_concordance(snp_labels, num_hets + 0.015, c_rare, 'r', 'Imputed Samples, Rare', threshold=threshold)

#    P.ylim([-0.05, 1.05])
    P.xlabel('# WGS Hets')
    P.ylabel('iPlex-Imputation Concordance')
    P.title('Imputation Accuracy for Rare Variants in the Hutterites')
    P.legend(loc='lower right', prop={'size': 10})

def hist_rare_imputation_concordance(c_rare, fig_num):
    '''Histogram of the rare samples imputation concordance.'''
    fig = P.figure(fig_num)
    P.clf() 
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    ax.hist(c_rare, bins=50, normed=True)
    ax.set_xlim([0, 1])
    ax.set_xlabel('iPlex-Imputation Concordance C')
    ax.set_ylabel('Frequency (%)')
    P.title('Imputation Accuracy for Samples with the Rare Allele')
    s = 'C = 1.0               : %5.1f%%\n' \
        '0.95 <= C < 1.0: %5.1f%%' % \
        (100.0 * len(P.mlab.find(c_rare == 1)) / len(c_rare),
         100.0 * len(P.mlab.find((c_rare >= 0.95) & (c_rare < 1))) / len(c_rare))
    ax.text(0.03, P.ylim()[1] * 0.97 - 3.3, s,
        bbox={'facecolor':'yellow', 'alpha':0.5, 'pad': 10})
 
def plot_concordance(labels, xdata, concordance, color, label, threshold=0.0):
    P.subplots_adjust(bottom=0.1)
    P.scatter(xdata, concordance, marker='^', s=100, lw=0, color=color, label=label)
    # c=concordance, cmap=P.get_cmap('Spectral'))
    # Annotate SNPs with low conocordance Annotate points with poor a
    for label, x, y in ((label, x, y) for (label, x, y) in zip(labels, xdata, concordance) if y < threshold):
        P.annotate(label, xy=(x, y), xytext=(80, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
