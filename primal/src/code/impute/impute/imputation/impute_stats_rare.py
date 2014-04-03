#!/usr/bin/env python
'''
============================================================
Calculate imputation statistics for the rare allele
study. Generate summary plots. TODO: merge with
impute_stats.

Created on December 3, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib.pyplot as P, numpy as np, util, impute as im, collections
from impute.tools import recode
from impute.data.constants import SMALL_FLOAT

#---------------------------------------------
# Methods - Statistics
#---------------------------------------------

'''Commonly-used filters to compare a pair of genotypes.'''
# Both genotypes are called
CALLED = lambda r1, r2: (r1 > 0) & (r2 > 0)

# At least one genotype contains the minor allele (2); the other might not be called
HAS_MINOR_ALLELE = lambda r1, r2: (r1 == 3) | (r1 == 4) | (r2 == 3) | (r2 == 4)

# Both genotypes are called and at least one contains the minor allele (2)
CALLED_AND_HAS_MINOR_ALLELE = lambda r1, r2: CALLED(r1, r2) & HAS_MINOR_ALLELE(r1, r2)


#---------------------
#class ImputationStats(object):
#    def __init__(t, snp=None):
#        self.t = t

def stats(t, snp=None):
    '''Return a record array with imputation statistics.'''
    T = t.sample_index_to_impute
    imputed = t.imputed_data[:, T, :]
    tot_to_impute = 2 * imputed.shape[1]
    snp = snp if snp is not None else np.arange(t.num_snps)
    stats = np.zeros((len(snp),),
                     dtype=[
                            ('dist_cm', 'f4'), # Genetic distance from beginning of chromosome
                            ('count', '(2,)i4'), # Allele count
                            ('frequency', '(2,)f4'), # Allele frequency
                            ('call_rate', 'f4'), # Imputation Call rate
                            ('call_rate_training', 'f4')  # Imputation Call rate
                            ])
    call_rate_training = 1.0 * np.sum(np.sum(t.imputed_data[:, t.sample_index, :] != 0, axis=2), axis=1)# / (2 * len(t.sample_index))        
    for row, snp_index in enumerate(snp):
        # TODO: replace by a bulk group-by/hist?
        # g = t.__t.training_data[snp_index, :, :]
        i = imputed[snp_index, :]
        (c1, c2) = (len(np.where(i == 1)[0]), len(np.where(i == 2)[0]))
        c = c1 + c2 + SMALL_FLOAT
        f1, f2 = (1.0 * c1) / c, (1.0 * c2) / c
        call_rate = 1.0 * len(i.nonzero()[0]) / tot_to_impute
        # print 'c1 %4d c2 %4d f1 %.2f f2 %.2f call rate %5.2f' % (c1, c2, f1, f2, call_rate)
        stats[row] = (t.snp['dist_cm'][snp_index], [c1, c2], [f1, f2], call_rate, call_rate_training[snp_index])
    return stats

def concordance(g1, g2, data_filter=CALLED, samples=None):
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
    maf = 100.0 * im.gt.allele_frequencies_by_snp(g)[1]

    c_training = concordance(iPlex, cgi, data_filter=CALLED)
    c_training_rare = concordance(iPlex, cgi, data_filter=HAS_MINOR_ALLELE)
    
    c_imputed = concordance(g, i, data_filter=CALLED, samples=imputed.sample_index_to_impute)
    c_rare = concordance(g, i, data_filter=CALLED_AND_HAS_MINOR_ALLELE, samples=imputed.sample_index_to_impute)
    
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

def plot_call_rate2(stats):
    # Run after run_impute_chr. TODO: consolidate with the stats object above
    x = stats['dist_cm']
    maf = np.min(stats['frequency'], axis=1)
    call_rate = stats['call_rate']
    call_rate_training = stats['call_rate_training']
    P.clf()
    P.hold(True)
    P.plot(x, maf, 'bo-', label='MAF')
    P.plot(x, call_rate, 'ro-', label='Call Rate, All')
    P.plot(x, call_rate_training, 'ro-', label='Call Rate, Training')
    
    P.legend(loc='upper left')
    P.xlabel('Genetic Distance [cM]')
    P.title('Chr 22 Imputation Results')
