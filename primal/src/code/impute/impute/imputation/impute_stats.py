#!/usr/bin/env python
'''
============================================================
Calculate imputation statistics; generate summary plots.

Created on December 3, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib.pyplot as P, numpy as np, util, impute as im, sys, db_gene
from impute.tools import recode
from impute.data.constants import SMALL_FLOAT, ALLELES, SNP, SAMPLE, NUM_CHROMOSOMES
from impute.plot.plots import genome_cm_offset, genome_bp_offset
from scipy import ndimage

#---------------------------------------------
# Constants
#---------------------------------------------
'''
Commonly-used filters to compare a pair of genotypes.
Recoding a genotype to a single integer: 
(0,0) ->  0
(0,1) -> -1
(0,2) -> -2
(1,1) ->  2
(1,2) ->  3
(2,2) ->  4
'''

# All SNPs
ALL = lambda r1, _: (r1 == r1)

# Both genotypes are fully called
CALLED = lambda r1, r2: (r1 > 0) & (r2 > 0)

# True genotypes are fully called, imputed are partially or fully called 
GENOTYPE_CALLED = lambda r1, r2: (r1 > 0) & (r2 != 0)

# Imputed is fully called
IMPUTED_CALLED = lambda r1, r2: (r2 > 0)

# At least one genotype contains the minor allele (2); the other might not be called
HAS_MINOR_ALLELE = lambda r1, r2: (r1 == 3) | (r1 == 4) | (r2 == 3) | (r2 == 4)

# Both genotypes are called and at least one contains the minor allele (2)
CALLED_AND_HAS_MINOR_ALLELE = lambda r1, r2: CALLED(r1, r2) & HAS_MINOR_ALLELE(r1, r2)

####################################################################################
class StatsComputer(object):
    '''Computes concordances between imputed results and a golden standard (known genotypes).'''
    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, t, g):
        '''t = ImputationSet imputed data set.
        g = Genotype data for all samples; sample indices must correspond to t''s sample indices.
        '''
        self.t = t
        self.T = t.sample_index  # training set
        self.I = t.sample_index_to_impute  # imputed sample set
        self.imputed_data = t.imputed_data[:, self.I, :]
        self.genetic_coord = self.t.snp['dist_cm']

        self.g = g
        self.rg = recode.recode_single_genotype(g)
        self.ri = recode.recode_single_genotype(t.imputed_data)
        self.maf = np.amin(im.gt.allele_frequencies_by_snp(g), axis=0)
        # self.concordance_training = self.concordance(SNP, samples=self.T)
        # self.concordance_imputed = self.concordance(SNP, samples=self.I)
        # self.concordance_all = self.concordance(SNP)
       
        # Lazily-initialized cached properties
        self._allele_count = None
        self._frequency = None
       
    #---------------------------------------------
    # Properties
    #---------------------------------------------    
    @staticmethod
    def _intersect_samples(x, y) :
        return np.intersect1d(x, y) if y is not None else x
    
    @property
    def allele_count(self):
        '''Return the allele count in the genotypes of all imputed samples that have been (partially
        or fully) called, broken down by SNP.'''
        if self._allele_count is None:
            self._allele_count = [np.sum(np.sum(self.t.imputed_data[:, self.I, :] == allele, axis=2), axis=1) for allele in ALLELES]        
        return self._allele_count
     
    @property
    def frequency(self):
        '''Return the allele frequency in the genotypes of all imputed samples that have been called,
        broken down by SNP.'''
        if self._frequency is None:
            count = self.allele_count
            tot = count[0] + count[1] + SMALL_FLOAT
            self._frequency = [(1.0 * c) / tot for c in count]
        return self._frequency

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def dist(self, s):
        '''Distance to lowest common ancestor between a sample and the training set.'''
        return min(im.pt.lowest_common_ancestor(self.t.pedigree.graph, s, x)[1] for x in self.t.sample_index)

    def sample_call_rate_and_concordance(self, snps):
        '''Return a tuple of call rate and concordance of samples at selected SNPs.''' 
        return np.array((self.call_rate_of_sample(snps=snps), self.concordance(SAMPLE, snps=snps))).transpose()

    def sample_concordance_array(self, sample, data_filter=CALLED, snps=None):
        '''Return a binary array indicating whether the sample is concordant or discordant at each
        SNP in the SNP list snps (or all SNPs, if snps=None).'''
        r1, r2, _, _ = self._restrict_arrays(data_filter, [sample], snps)
        return r1 == r2

    def worst_sample_concordance(self, n, data_filter=CALLED, snps=None):
        '''Return the concordance by SNP among the most discordant n samples.'''
        c = self.concordance(SAMPLE, data_filter=data_filter, snps=snps)
        worst = np.argsort(c)[0:n]
        return worst, self.concordance(SNP, data_filter=data_filter, samples=worst, snps=snps)
    
    def call_rate_imputed(self, data_filter=CALLED, snps=None, samples=None, counts=False):
        '''Return the call rate of each SNP over all imputed samples.
        If counts=True, returns counts instead of call rates.'''
        functor = self._call_counts_by_axis if counts else self._call_rate 
        return functor(SNP, data_filter=data_filter, samples=StatsComputer._intersect_samples(self.I, samples), snps=snps)

    def call_rate_training(self, data_filter=CALLED, snps=None, samples=None, counts=False):
        '''Return the call rate (= imputation-phasing %) of each SNP in the training sample set.
         If counts=True, returns counts instead of call rates.'''
        functor = self._call_counts_by_axis if counts else self._call_rate 
        return functor(SNP, data_filter=data_filter, samples=StatsComputer._intersect_samples(self.T, samples), snps=snps)
   
    def call_rate_of_sample(self, data_filter=CALLED, snps=None, samples=None, counts=False):
        '''Return the call rate of each sample over SNPs i.
        If counts=True, returns counts instead of call rates.'''
        functor = self._call_counts_by_axis if counts else self._call_rate
        return functor(SAMPLE, data_filter=data_filter, snps=snps, samples=samples)

    def genotype_call_rate_of_sample(self, data_filter=CALLED, snps=None, samples=None):
        '''Return the call rate of each sample over SNPs i.
        If counts=True, returns counts instead of call rates.'''
        functor = self._genotype_call_rate
        return functor(SAMPLE, data_filter=data_filter, snps=snps, samples=samples)

    def concordance(self, axis, data_filter=CALLED, samples=None, snps=None, counts=False):
        '''Return the concordance (# of filled entries that agree) between two genotype data arrays g1,g2.
        Returns an array of SNP indices that have data for comparison, and the corresponding concordance
        rates. If counts=True, returns counts instead of call rates.'''
        r1_called, r2_called, called, r1_shape = self._restrict_arrays(data_filter, samples, snps)
        return StatsComputer._concordance_counts_by_axis(r1_called, r2_called, called, r1_shape[axis], axis)[1:3] if counts \
            else StatsComputer._concordance_by_axis(r1_called, r2_called, called, r1_shape[axis], axis)

    def compare_at_snp(self, snp):
        '''Compare Imputed with original genotypes at snp # snp.''' 
        return np.concatenate((np.arange(self.t.num_samples)[np.newaxis].transpose(),
                               self.t.imputed_data[snp], self.g[snp]), axis=1)

    def stats_snp(self, snps=None, samples=None, counts=False):
        '''Calculate all statistics vs. SNP of interest without plotting.'''
        snps = snps if snps is not None else np.arange(self.t.num_snps)
        x = self.genetic_coord[snps]
        call_rate_imputed_full = self.call_rate_imputed(snps=snps, samples=samples, counts=counts)
        call_rate_imputed_partial = self.call_rate_imputed(snps=snps, samples=samples , data_filter=GENOTYPE_CALLED, counts=counts)
        call_rate_training = self.call_rate_training(snps=snps, samples=samples, counts=counts)
        imputed_samples = StatsComputer._intersect_samples(self.I, samples)
        concordance_imputed = self.concordance(SNP, samples=imputed_samples, counts=counts)[snps]
        concordance_imputed_het = self.concordance(SNP, samples=imputed_samples, counts=counts, data_filter=CALLED_AND_HAS_MINOR_ALLELE)[snps]
        return snps, x, call_rate_imputed_full, call_rate_imputed_partial, call_rate_training, concordance_imputed, concordance_imputed_het

    def stats_maf(self, snps, call_rate_imputed_full, call_rate_training, concordance_imputed, concordance_imputed_het):
        '''Calculate all statistics vs. MAF of interest without plotting.'''
        maf = self.maf[snps]
        bins = np.arange(min(maf), max(maf), 0.02)
        r = range(1, len(bins))
        digitize = lambda x: np.digitize(x, bins)
        maf_bin = digitize(maf)
        avg = lambda x: np.array(map(lambda m: np.mean(x[np.where(maf_bin == m)[0]]), r))
        call_rate_imputed = avg(call_rate_imputed_full)
        has_data = np.where(np.isfinite(call_rate_imputed))[0]
        call_rate_imputed = call_rate_imputed[has_data]
        maf = avg(maf)[has_data]
        call_rate_training = avg(call_rate_training)[has_data]
        concordance_imputed = avg(concordance_imputed)[has_data]
        concordance_imputed_het = avg(concordance_imputed_het)[has_data]
        return maf, call_rate_imputed, call_rate_training, concordance_imputed, concordance_imputed_het
        
    def stats_sample(self, snps=None, samples=None, counts=False):
        '''Calculate all statistics vs. sample of interest without plotting.'''
        call_rate_full = self.call_rate_of_sample(snps=snps, samples=samples, counts=counts)
        # Allele call rate = (call_rate_full + call_rate_partial)/2
        call_rate_partial = self.call_rate_of_sample(snps=snps, data_filter=GENOTYPE_CALLED, samples=samples, counts=counts)
        call_rate_partial = 0.5 * (call_rate_partial + call_rate_full)
        concordance = np.array(self.concordance(SAMPLE, snps=snps, samples=samples, counts=counts))
        return call_rate_full, call_rate_partial, concordance
        
    def stats(self, snps=None, samples=None, counts=False):
        '''Calculate all statistics of interest in a Stats object.'''
        result = Stats()
        result.snp.info = self.t.snp
        result.sample.samples = samples if samples is not None else np.arange(self.g.shape[1])   
        result.snp.snps, result.snp.x, result.snp.call_rate_imputed_full, \
        result.snp.call_rate_imputed_partial, result.snp.call_rate_training, \
        result.snp.concordance_imputed, result.snp.concordance_imputed_het = \
        self.stats_snp(snps=snps, samples=samples, counts=counts)
        result.maf.maf, result.maf.call_rate_imputed, result.maf.call_rate_training, \
        result.maf.concordance_imputed, result.maf.concordance_imputed_het = \
        self.stats_maf(result.snp.snps, result.snp.call_rate_imputed_full, result.snp.call_rate_training, result.snp.concordance_imputed, result.snp.concordance_imputed_het)
        result.sample.call_rate_full, result.sample.call_rate_partial, result.sample.concordance = self.stats_sample(snps=snps, samples=samples)
        result.snp.maf = self.maf[snps]
        return result

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def _restrict_arrays(self, data_filter=CALLED, samples=None, snps=None):
        '''Return an array of SNP indices that have data for concordance comparison between self.rg, self.ri.'''
        r1, r2 = self.rg, self.ri
        if samples is not None:
            if len(samples) == 0: r1, r2 = np.zeros((r1.shape[0], 0)), np.zeros((r2.shape[0], 0))
            else: r1, r2 = r1[:, samples], r2[:, samples]
        if snps is not None and r1.size: r1, r2 = r1[snps, :], r2[snps, :]
        called = np.where(data_filter(r1, r2))
        return r1[called], r2[called], called, r1.shape
    
    def _call_rate(self, axis, data_filter=CALLED, samples=None, snps=None):
        '''Return the imputation call rate by axis (0=SNP, 1=sample).'''
        _, ri, called, num_groups = self._restrict_arrays(data_filter, samples, snps)
        num_called, has_data = util.sum_by_group(np.ones_like(ri).astype(np.byte), called[axis])
        call_rate = np.zeros((num_groups[axis],), dtype=np.float)
        call_rate[has_data] = (1.0 * num_called) / num_groups[1 - axis]
        return call_rate

    def _genotype_call_rate(self, axis, data_filter=CALLED, samples=None, snps=None):
        '''Return the call rate by axis (0=SNP, 1=sample).'''
        rg, _, called, num_groups = self._restrict_arrays(data_filter, samples, snps)
        num_called, has_data = util.sum_by_group(np.ones_like(rg).astype(np.byte), called[axis])
        call_rate = np.zeros((num_groups[axis],), dtype=np.float)
        call_rate[has_data] = (1.0 * num_called) / num_groups[1 - axis]
        return call_rate

    def _call_counts_by_axis(self, axis, data_filter=CALLED, samples=None, snps=None):
        '''Return the call counts (called count, all count) by axis (0=SNP, 1=sample).'''
        _, ri, called, num_groups = self._restrict_arrays(data_filter, samples, snps)
        num_called, has_data = util.sum_by_group(np.ones_like(ri).astype(np.byte), called[axis])
        call_rate = np.zeros((num_groups[axis],), dtype=np.float)
        num_all = num_groups[1 - axis]
        call_rate[has_data] = (1.0 * num_called) / num_all
        # return has_data, num_called, num_all
        return num_called, num_all

    @staticmethod
    def _concordance_by_axis(r1, r2, groups, num_groups, axis):
        '''Calculate concordances by axis (0=SNP, 1=sample).'''
        # print np.concatenate((g1[groups], g2[groups]), axis=1)
        has_data, _, _, concordant, discordant = StatsComputer._concordance_counts_by_axis(r1, r2, groups, num_groups, axis)
        c = np.ones((num_groups,), dtype=np.float)  # no data => concordance=1 by convention
        c[has_data] = (1.0 * concordant) / (concordant + discordant)
        return c

    @staticmethod
    def _concordance_counts_by_axis(r1, r2, groups, num_groups, axis):
        '''Calculate #concordances, #discordances by axis (0=SNP, 1=sample).'''
        # print np.concatenate((g1[groups], g2[groups]), axis=1)
        groups = groups[axis]
        concordant, has_data = util.sum_by_group((r1 == r2).astype(np.byte), groups)
        discordant, _ = util.sum_by_group((r1 != r2).astype(np.byte), groups)
        c = np.zeros((num_groups,), dtype=np.float)
        c[has_data] = concordant
        d = np.zeros((num_groups,), dtype=np.float)
        d[has_data] = discordant
        return has_data, c, d, concordant, discordant
    
####################################################################################
class Stats(object):
    '''Holds statistical data of the comparison imputed results and a golden standard (known genotypes).
    This is a stand-alone data structure from which all plots can be generated without access to
    the original data This object is produced by StatsComputer.'''
    class _Struct(object): pass
    
    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self):
        '''Initialize an empty stats holder.'''
        self.snp = Stats._Struct()
        self.maf = Stats._Struct()
        self.sample = Stats._Struct()
        self._chrom = None
        self.linewidth = 1
    
    #---------------------------------------------
    # Properties - Reporting
    #---------------------------------------------
    @property
    def chrom(self):
        if self._chrom is None:
            snp_chrom = self.snp.info['chrom']
            self._chrom = snp_chrom[0] if len(np.unique(snp_chrom)) == 1 else 0
        return self._chrom
        
    @property
    def cm_cumulative(self):
        '''cM axis offset for aligning multiple chromosomes.'''
        if self.chrom: return self.snp.x
        else:
            offset = genome_cm_offset()[:NUM_CHROMOSOMES]
            return self.snp.x + offset[self.snp.info['chrom'] - 1]

    @property
    def bp_cumulative(self):
        '''cM axis offset for aligning multiple chromosomes.'''
        offset = genome_bp_offset()[:NUM_CHROMOSOMES]
        return self.snp.info['base_pair'] + offset[self.snp.info['chrom'] - 1]
    
    @property
    def cm_edge_dist(self):
        '''Distance to nearest edge of chromosome [cM].'''
        l = np.array(db_gene.snp.file_dao.DEFAULT_FILE_DAOS.chrom_dao.total_cm())
        return np.minimum(self.snp.x, l[self.chrom - 1] - self.snp.x) if self.chrom else \
            np.minimum(self.snp.x, l[self.snp.info['chrom'] - 1] - self.snp.x)

    #---------------------------------------------
    # Methods - Reporting
    #---------------------------------------------
    def samples_with_call_rates_above(self, threshold):
        '''Return samples with high call rates.'''
        t = self.sample
        return t.samples[np.where(t.call_rate_full >= threshold)[0]]
        
    def summary(self, out=sys.stdout):
        '''Print a summary report to the output stream out.'''
        snps = self.snp.snps
        samples = self.sample.samples
        out.write('Validation Results, #snps = %d, #samples = %d\n' % (len(snps), len(samples)))

        labels = ['Call Rate, Imputed (SNP)',
                  'Partial Call Rate, Imputed (SNP)',
                  'Phasing Rate, WGS (SNP)',
                  'Concordance (SNP)',
                  'Het Concordance (SNP)',
                  'Call Rate (sample)',
                  'Partial Call Rate (sample)',
                  'Concordance (sample)']
        for k, x in enumerate([self.snp.call_rate_imputed_full,
                               self.snp.call_rate_imputed_partial,
                               self.snp.call_rate_training,
                               self.snp.concordance_imputed,
                               self.snp.concordance_imputed_het,
                               self.sample.call_rate_full,
                               self.sample.call_rate_partial,
                               self.sample.concordance]):
            out.write('%-32s  %.4f +- %.4f  median %.4f  min %.4f  max %.4f\n' % \
            (labels[k], np.mean(x), np.std(x), np.median(x), np.min(x), np.max(x)))
        worst_snp = np.argmin(self.snp.concordance_imputed)
        out.write('Most discordant SNP = %d at %.4f\n' % (snps[worst_snp], self.snp.concordance_imputed[worst_snp]))

    #---------------------------------------------
    # Methods - Plots
    #---------------------------------------------    
    def plot_vs_snp(self, snp_style='continuous', filter_length=20, x_axis='cm_cumulative'):
        '''Plot call rates and imputation accuracy vs. SNP.'''
        if x_axis == 'cm_cumulative':
            x = self.cm_cumulative
            xlabel = 'Genetic Distance from Beginning of Genome [cM]'
        elif x_axis == 'bp_cumulative':
            x = self.bp_cumulative
            xlabel = 'Physical Position from Beginning of Genome [bp]'
        elif x_axis == 'cm_edge_dist':
            sz = 2
            x = self.cm_edge_dist
            xlabel = 'Genetic Distance from Nearest Chromosome Edge [cM]'
        r = self.snp
        k = np.argsort(x)

        P.clf()
        P.hold(True)
        self._draw_chrom_dividers(x_axis)

        # Smoothing filter
        f = lambda d: ndimage.median_filter(d, filter_length)[k]
        
        if snp_style == 'continuous':
            x = x[k]
            if x_axis == 'cm_edge_dist':
                P.scatter(x, f(r.call_rate_imputed_full), lw=0, s=sz, color='g', label='Genotype Call, Imputed', linewidth=self.linewidth)
                P.scatter(x, f(r.call_rate_imputed_partial), lw=0, s=sz, color='b', label='Allele Call, Imputed', linewidth=self.linewidth)
                P.scatter(x, f(r.call_rate_training), lw=0, s=sz, color='c', label='Phasing Rate, Training', linewidth=self.linewidth)
                P.scatter(x, f(r.concordance_imputed), lw=0, s=sz, color='r', label='Concordance', linewidth=self.linewidth)
                P.scatter(x, f(r.concordance_imputed_het), lw=0, s=sz, color='m', label='Het Concordance', linewidth=self.linewidth)
            else:
                P.plot(x, f(r.call_rate_imputed_full), 'g-', label='Genotype Call, Imputed', linewidth=self.linewidth)
                P.plot(x, f(r.call_rate_imputed_partial), 'b-', label='Allele Call, Imputed', linewidth=self.linewidth)
                P.plot(x, f(r.call_rate_training), 'c-', label='Phasing Rate, Training', linewidth=self.linewidth)
                P.plot(x, f(r.concordance_imputed), 'r-', label='Concordance', linewidth=self.linewidth)
                P.plot(x, f(r.concordance_imputed_het), 'm-', label='Het Concordance', linewidth=self.linewidth)
        else:
            snp_labels = map(lambda x: 'chr%d:%s' % (x['chrom'], x['name'].split(':')[-1]), self.snp.info)
            self._scatter(snp_labels, x, r.call_rate_imputed_full, 'g', 'Genotype Call, Imputed', threshold=0.5)
            self._scatter(snp_labels, x, r.call_rate_imputed_partial, 'b', 'Allele Call, Imputed', threshold=0.5)
            self._scatter(snp_labels, x, r.call_rate_training, 'c', 'Phasing Rate, Training', threshold=0.5)
            self._scatter(snp_labels, x, r.concordance_imputed, 'r', 'Concordance', threshold=0.5)
            self._scatter(snp_labels, x, r.concordance_imputed_het, 'm', 'Het Concordance', threshold=0.5)
        
        P.xlim([min(x) - 5, max(x)])
        # P.ylim([0.9 * min(r.call_rate_imputed_full), 1.02])
        P.ylim([-0.02, 1.02])
        #P.legend(loc='best')
        P.legend(loc='lower right')
        P.xlabel(xlabel)
        P.ylabel('Call Rate, Concordance')
        P.title('Imputation Performance' + self._chrom_suffix() + '; filter size = ' + repr(filter_length))
            
    def scatter_snp_concordance(self, snp_style='continuous'):
        '''Scatter of SNP het concordance vs. concordance.'''
        c = self.snp.concordance_imputed
        ch = self.snp.concordance_imputed_het
        P.clf()
        P.hold(True)
        P.grid(True)
        P.scatter(c, ch, color='b', lw=0)
        P.xlim([min(c) - 0.02, 1.02])
        P.ylim([min(ch) - 0.02, 1.02])
        x = np.linspace(0, 1, 100)
        if snp_style == 'continuous': P.plot(x, x, 'k--')
        P.xlabel('Concordance')
        P.ylabel('Het Concordance')
        P.title('Het concordance vs. concordance, $R^2$ = %.4f' % 
                (np.corrcoef(c + SMALL_FLOAT, ch + SMALL_FLOAT)[0, 1] ** 2,))

    def plot_vs_maf(self, snp_style='continuous'):
        '''Plot call rates and imputation accuracy vs. minor allele frequency.'''
        P.clf()
        P.hold(True)
        ax = P.axes()
        if snp_style == 'continuous':
            r = self.maf
            ax.plot(r.maf, r.concordance_imputed, 'r-', label='Concordance (all genotypes)', linewidth=self.linewidth)
            ax.plot(r.maf, r.concordance_imputed_het, 'g-', label='Concordance (het genotypes)', linewidth=self.linewidth)
            ax.plot(r.maf, r.call_rate_imputed, 'b-', label='Call Rate, Imputed', linewidth=self.linewidth)
            if r.call_rate_imputed.size: P.ylim([0.99 * min(r.call_rate_imputed), 1.02])
        else:
            r = self.snp
            ax.plot(r.maf, r.concordance_imputed, 'r^', label='Concordance (all genotypes)', linewidth=self.linewidth)
            ax.plot(r.maf, r.concordance_imputed_het, 'g^', label='Concordance (het genotypes)', linewidth=self.linewidth)
            ax.plot(r.maf, r.call_rate_imputed_full, 'b^', label='Call Rate, Imputed', linewidth=self.linewidth)
#            snp_labels = map(lambda x: 'chr%d:%s' % (x['chrom'], x['name'].split(':')[-1]), self.snp.info)
#            self._scatter(snp_labels, r.maf, r.call_rate_imputed, 'g', 'Call Rate, Imputed', threshold=0.0)
#            self._scatter(snp_labels, r.maf, r.call_rate_training, 'b', 'Phasing Rate, Training', threshold=0.0)
#            self._scatter(snp_labels, r.maf, r.concordance_imputed, 'r', 'Concordance', threshold=0.0)

        # P.legend(loc='lower left', prop={'size': 20})
        P.legend(loc='lower right')
        P.xlabel('MAF')
        # P.ylabel('Call Rate, Concordance')
        # P.title('Imputation Performance' + self._chrom_suffix())
        # for item in ([ax.xaxis.label, ax.yaxis.label]):
        #   item.set_fontsize(30)
        # for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            # item.set_fontsize(20)

    def plot_vs_sample(self, snps=None, samples=None):
        '''Plot concordances of samples, sorted.'''
        r = self.sample
        i = np.argsort(r.call_rate_full)
                
        # ind = min(np.where(sorted_concordance >= 0.96)[0])
        P.clf()
        P.hold(True)
        P.plot(r.call_rate_full[i], 'g-', label='Call Rate', linewidth=self.linewidth)
        P.plot(r.call_rate_partial[i], 'b-', label='Partial Call Rate', linewidth=self.linewidth)
        P.plot(r.concordance[i], 'r-', label='Concordance', linewidth=self.linewidth)
        
        # P.xlim([-1, ind + 5])
        P.xlim([-20, len(r.concordance) + 1])
        ylim = P.ylim()
        P.ylim([ylim[0], 1.02])
        # P.plot(x, concordance_training, 'ro-', label='Training')
        # P.legend(loc='upper right')
        P.xlabel('Sample #')
        P.ylabel('Call Rate / Concordance')
        P.title('Imputation Performance' + self._chrom_suffix())
        P.legend(loc='lower right')
        # i = np.argsort(concordance)
        # return call_rate_full, call_rate_partial, concordance, np.array([(x, call_rate_full[x], call_rate_partial[x], concordance[x]) for x in i])
        # Add distance to nearest relative. Slow. 
#        return call_rate_full, call_rate_partial, concordance, \
#            np.array([(x, self.dist(x), call_rate_full[x], call_rate_partial[x], concordance[x])
#                      for x in i])

    @staticmethod
    def load(npz_file):
        '''Deserialize to object from an npz file.'''
        files = np.load(npz_file)
        stats = Stats()
        
        stats._chrom = int(files['chrom'])
        
        stats.sample.call_rate_full = files['sample_call_rate_full']
        stats.sample.concordance = files['sample_concordance']
        stats.sample.call_rate_partial = files['sample_call_rate_partial']  
        stats.sample.samples = files['sample_samples']
        
        stats.snp.call_rate_imputed_full = files['snp_call_rate_imputed_full']
        stats.snp.concordance_imputed_het = files['snp_concordance_imputed_het']
        stats.snp.call_rate_imputed_partial = files['snp_call_rate_imputed_partial']
        stats.snp.info = files['snp_info']
        stats.snp.call_rate_training = files['snp_call_rate_training']         
        stats.snp.snps = files['snp_snps']
        stats.snp.concordance_imputed = files['snp_concordance_imputed']        
        stats.snp.x = files['snp_x']
        stats.snp.maf = files['snp_maf'][0]
                  
        stats.maf.call_rate_imputed = files['maf_call_rate_imputed']        
        stats.maf.concordance_imputed_het = files['maf_concordance_imputed_het']
        stats.maf.call_rate_training = files['maf_call_rate_training'] 
        stats.maf.maf = files['maf_maf']
        stats.maf.concordance_imputed = files['maf_concordance_imputed']
                          
        return stats
        
    def save(self, npz_file):
        '''Serialize to object to an npz file.'''
        np.savez(npz_file,
                 chrom=self.chrom,
                 sample_call_rate_full=self.sample.call_rate_full,
                 sample_concordance=self.sample.concordance,
                 sample_call_rate_partial=self.sample.call_rate_partial,
                 sample_samples=self.sample.samples,
                 
                 snp_call_rate_imputed_full=self.snp.call_rate_imputed_full,
                 snp_concordance_imputed_het=self.snp.concordance_imputed_het,
                 snp_call_rate_imputed_partial=self.snp.call_rate_imputed_partial,
                 snp_info=self.snp.info,
                 snp_call_rate_training=self.snp.call_rate_training,
                 snp_snps=self.snp.snps,
                 snp_concordance_imputed=self.snp.concordance_imputed,
                 snp_x=self.snp.x,
                 snp_maf=self.snp.maf,
                  
                 maf_call_rate_imputed=self.maf.call_rate_imputed,
                 maf_concordance_imputed_het=self.maf.concordance_imputed_het,
                 maf_call_rate_training=self.maf.call_rate_training,
                 maf_maf=self.maf.maf,
                 maf_concordance_imputed=self.maf.concordance_imputed,
                 # sample_index=self.sample_index,
                 # pedigree=np.array([self.pedigree])
                 )

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------    
    def _scatter(self, labels, xdata, ydata, color, plot_label, threshold=0.0):
        P.subplots_adjust(bottom=0.1)
        P.scatter(xdata, ydata, marker='^', s=100, lw=0, color=color, label=plot_label)
        # c=concordance, cmap=P.get_cmap('Spectral'))
        # Annotate SNPs with low conocordance Annotate points with poor a
        for label, x, y in ((label, x, y) for (label, x, y) in zip(labels, xdata, ydata) if y < threshold):
            P.annotate(label, xy=(x, y), xytext=(80, 20),
                textcoords='offset points', ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    def _draw_chrom_dividers(self, x_axis):
        if self.chrom == 0:
            if x_axis == 'cm_cumulative':
                genome_offset = genome_cm_offset()
            else:
                genome_offset = genome_bp_offset()
            offset = genome_offset[:NUM_CHROMOSOMES]
            for chr_offset in offset[1:]: P.plot([chr_offset, chr_offset], [-0.1, 1.1], c='k', linestyle=':')

    def _chrom_suffix(self):
        return ((': Chromosome %d' % (self.chrom,)) if self.chrom else '')

####################################################################################
#---------------------------------------------
# Methods - CGI Imputation Stats
#---------------------------------------------

def read_call_rates(count_file, num_samples):
    count_called = np.loadtxt(count_file, usecols=range(1, 27, 3), dtype=np.long)

    chrom = range(1, count_called.shape[0] + 1)  # Chromosome number
    count_full = np.sum(count_called[:, np.array([0, 1, 3, 4])], axis=1)
    count_all = np.sum(count_called, axis=1)
    count_partial = count_all - count_called[:, 8]
    
    call_rate_full = (1.0 * count_full) / count_all
    call_rate_partial = (1.0 * count_partial) / count_all

    total_full = np.sum(count_full, axis=0)
    total_partial = np.sum(count_partial, axis=0)
    total_all = np.sum(count_all)
    
    return util.Struct(chrom=chrom,
                       call_rate_full=call_rate_full,
                       call_rate_partial=call_rate_partial,
                       num_snps=(1.0 * total_all) / num_samples, \
                       total_call_rate=(100. * total_full) / total_all, \
                       total_call_rate_partial=(100. * total_partial) / total_all)

def plot_cgi_call_rate(count_file, count_phasing_file, num_samples, title=None):
    '''Plot CGI call rates.'''
    samples = read_call_rates(count_file, num_samples)
    phasing = read_call_rates(count_phasing_file, num_samples)

    P.figure(1)
    P.clf()
    P.hold(True)
    P.plot(samples.chrom, samples.call_rate_full, 'ro-', label='Genotype')
    P.plot(samples.chrom, samples.call_rate_partial, 'bo-', label='Allele')
    P.plot(samples.chrom, phasing.call_rate_full, 'go-', label='Phasing')
    P.legend(loc='upper right', prop={'size': 10})
    P.xlabel('Chromosome #')
    P.ylabel('Call Rate')
    P.xlim([0, 23])
    P.ylim([0.4, 1])
    title = title if title else 'Hutterites CGI Imputation Call Rates'
    P.title('%s\nN = %d, Imputed loci = %d\nOverall Call Rate: Genotype %.1f%% Allele %.1f%%' % \
            (title, num_samples, samples.num_snps, samples.total_call_rate, samples.total_call_rate_partial))

    print 'N = %d, Imputed loci = %d, Phasing %.1f%% Genotype %.1f%% Allele %.1f%%' % \
    (num_samples, samples.num_snps, phasing.total_call_rate, samples.total_call_rate, samples.total_call_rate_partial)
    
    P.show()
