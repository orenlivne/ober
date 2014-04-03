#!/usr/bin/env python
'''
============================================================
Plot the results of a phasing batch pipeline. 
 
Created on September 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, numpy as np, matplotlib.pyplot as P, impute as im, os

class BatchStats(object):
    def __init__(self, stats):
        s = np.concatenate([x for x in stats])
        self.s = s
        sum_all_parts = lambda data: reduce(lambda x, y: x + y, (data(x) for x in s))
        self.num_samples = len(sum_all_parts(lambda x: x.imputed.by_sample))
        self.total_genotypes = (1.0 * sum(x.num_genotypes for x in s)) / self.num_samples
        self.samples = np.arange(1, self.num_samples + 1)

        self.missing_orig = sum_all_parts(lambda x: x.num_snps - x.called_orig.by_sample)
        self.imputed = sum_all_parts(lambda x: x.imputed.by_sample)
        self.imputed_partial = sum_all_parts(lambda x: x.imputed_partial.by_sample)
        self.called_orig_fraction = sum_all_parts(lambda x: x.called_orig.by_sample) / self.total_genotypes
        self.called_fraction = sum_all_parts(lambda x: x.called_orig.by_sample + x.imputed.by_sample) / self.total_genotypes
        
        self.missing_orig_fraction = 100.0 * self.missing_orig / self.total_genotypes
        self.imputed_fraction = (100.0 * self.imputed) / self.missing_orig
        self.imputed_partial_fraction = (100.0 * self.imputed_partial) / self.missing_orig
        
        self.imputed_by_chr = np.array([sum(x.imputed.total for x in c) for c in stats])
        self.errors_by_chr = np.array([sum(x.errors.total for x in c) for c in stats])
        self.called_orig_by_chr = np.array([sum(x.called_orig.total for x in c) for c in stats])
        self.total_by_chr = np.array([sum(x.num_genotypes for x in c) for c in stats])
        self.missing_by_chr = self.total_by_chr - self.called_orig_by_chr - self.imputed_by_chr - self.errors_by_chr
        
        self.imputed_fraction_by_chr = (100.0 * self.imputed_by_chr) / self.total_by_chr
        self.called_orig_fraction_by_chr = (100.0 * self.called_orig_by_chr) / self.total_by_chr

def plot_imputation_results(s):
    '''A big plot of imputation results - summary statistics.'''
    P.suptitle('Imputation Results, N=%d' % (s.num_samples,), fontsize=16)
    
    P.subplot(221)
    P.plot(s.num_samples - s.samples + 1, np.sort(s.called_orig_fraction), 'b', label='Original')
    P.hold(True)
    P.plot(s.num_samples - s.samples + 1, np.sort(s.called_fraction), 'r', label='Imputed')
    P.xlim([900, s.num_samples])
    P.grid()
    P.xlabel('# Retained Samples')
    P.ylabel('Min. Call Rate (Fill %)')
    P.title('Call Rate by Sample (Sorted Ascending)')
    P.legend(loc='upper right', prop={'size': 10})

    P.subplot(222)
    plot_missing_by_platform(s.missing_orig_fraction, s.imputed_fraction, s.pedigree, s.platform, s.colors)
#    P.scatter(s.missing_orig_fraction, s.imputed_fraction)
#    P.xlabel('% Originally Missing Genotypes')
#    P.ylabel('% Imputed Genotypes')
#    P.title('% Imputed vs. % Originally Missing by Sample')
#    P.ylim([-5, 105])
#    (_, xmax) = P.xlim()
#    P.xlim([-1, xmax])
    
#    P.hist(s.imputed_fraction, 30, normed=True, histtype='bar', color=['blue'], label=['% Imputed Genotypes'])
#    P.title('Imputed Genotype Percentage by Sample')
#    P.ylabel('Frequency in Samples')
#    P.xlabel('# Imputed / # Missing')
#    P.xlim((0, 100))
#    P.legend(loc='upper left', prop={'size': 10})

    P.subplot(223)
    P.hist([s.imputed_partial_fraction, s.imputed_fraction], 10, normed=True, histtype='bar',
           color=['red', 'blue'], label=['Fully Imputed', 'Partially Imputed'])    
    P.title('% Imputed Genotype by Sample')
    P.ylabel('Frequency in Samples')
    P.xlabel('% Imputed Genotypes')
    P.xlim((0, 100))
    P.legend(loc='upper right', prop={'size': 10})
    
    P.subplot(224)
    N = len(s.called_orig_fraction_by_chr)
    ind = np.arange(1, N + 1)  # the x locations for the groups
    width = 0.5  # the width of the bars: can also be len(x) sequence
    p1 = P.bar(ind, s.called_orig_fraction_by_chr, width, color='r')
    p2 = P.bar(ind, s.imputed_fraction_by_chr, width, color='y', bottom=s.called_orig_fraction_by_chr)
    P.title('% Imputed Genotypes by Chromosome')
    P.ylabel('% Genotypes')
    P.xlabel('Chromosome')
    P.xticks(ind[1::2] - 3 * width / 2., np.arange(1, N + 1, 2))
    # P.yticks(np.arange(0, 81, 10))
    P.legend((p1[0], p2[0]), ('% Originally Called', '% Imputed'), loc='lower right', prop={'size': 10})
    P.ylim([np.max((0, 100 - 3 * min(s.imputed_fraction_by_chr))), 100])

def hist_missing(missing_orig_fraction):
    '''Plot a histogram of missing genotypes in dthe data set.'''
    P.hist(s.missing_orig_fraction, 30, histtype='bar', color=['blue'], label=['% Missing Genotypes'])    
    P.xlabel('% Missing Genotypes')
    P.ylabel('Frequency')
    P.title('%% Originally Missing Genotypes (N = %d)' % (s.num_samples,))
    
def plot_missing_by_platform(missing_orig_fraction, imputed_fraction, pedigree, platform, colors):
    '''Compare missing data in different platforms.
    Per discussion with Dan, Carole, Mark on missing data distribution - why is there a gap
    between 2% and 6% in hist_missing()? Owing to different platforms? Answer: platform differences.'''
    c = np.zeros((pedigree.num_genotyped, 3), dtype=np.uint8)
    for (k, v) in platform.iteritems():
        ind = np.array([pedigree.node_of[x] for x in v if pedigree.is_genotyped(pedigree.node_of[x])])
        c[ind, :] = np.tile(colors[k], (len(ind), 1))

    P.scatter(s.missing_orig_fraction, s.imputed_fraction, c=c, norm=None)
    P.xlabel('% Originally Missing Genotypes')
    P.ylabel('% Imputed Genotypes')
    P.title('% Imputed vs. % Originally Missing by Sample and Platform')
    P.ylim([-5, 105])
    (_, xmax) = P.xlim()
    P.xlim([-0.2, xmax])
    P.legend([P.Rectangle((0, 0), 1, 1, fc=colors[k]) for k in colors.keys()], colors.keys(),
             loc='upper right', prop={'size': 10})

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # old_printoptions = np.get_printoptions()
    
    # Load data
    stats = np.load(sys.argv[1])['stats'][0]
    s = BatchStats(stats)
    # Compare platforms
    s.pedigree = im.io_pedigree.read(os.environ['OBER'] + '/data/lung/hutt-lung.pdg.tfam', genotyped_id_file=os.environ['OBER'] + '/data/lung/out-lung-1106/split_chr/chr1/hutt-lung_chr1_part0.tfam')
    e = np.loadtxt(os.environ['OBER'] + '/doc/findiv/excludeFindiv')
    s.platform = {
            '6.0': np.setdiff1d(np.loadtxt(os.environ['OBER'] + '/doc/findiv/6.0Findiv'), e),
            '5.0': np.setdiff1d(np.loadtxt(os.environ['OBER'] + '/doc/findiv/5.0Findiv'), e),
            '500k': np.setdiff1d(np.loadtxt(os.environ['OBER'] + '/doc/findiv/500Findiv'), e)
            }
    s.colors = {'6.0': [1, 0, 0], '5.0': [0, 1, 0], '500k': [0, 0, 1]}

    # Main imputation result
    P.close(1)
    P.figure(1, figsize=(13, 10), dpi=80)
    P.clf()
    plot_imputation_results(s)
    P.show()
    # P.savefig('imputation_results.png')
    P.savefig(sys.argv[2])
    
#    P.figure(2)
#    P.clf()
#    hist_missing(s.missing_orig_fraction)
#    P.savefig('missing.png')
#    P.figure(2)
#    P.clf()
#    plot_missing_by_platform(s.missing_orig_fraction, s.imputed_fraction, pedigree, platform, colors)
#    P.show()
#    P.savefig('missing_by_platform.png')
    
    # util.set_printoptions(old_printoptions)
