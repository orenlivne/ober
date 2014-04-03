#!/usr/bin/env python
'''
============================================================
Determine parent-of-origin for quasi-founders (genotyped
samples whose parents are not genotyped). The algorith,
compares the kinship of the father and the mother to each
of the IBD cliques containing the proband's hapotypes. 

Created on July 9, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, os, matplotlib.pyplot as P
from scipy import polyval
from stats.robust import linefit
from numpy.core.function_base import linspace
from scipy.stats.mstats_basic import mquantiles
from statutil import POSITIVE_TO_MINUS_1_TO_1
# from multiprocessing import Manager, Pool

#---------------------------------------------
# Constants
#---------------------------------------------
# Scaling factor for separation angle
ANGLE_SCALING_FACTOR = 2.0 / np.pi

#---------------------------------------------
# Methods
#---------------------------------------------
def plot_flip_measure_vs_sample(chrom, m):
    P.figure(2)
    P.clf()
    #P.title('Separation Measure M(C), Sorted: Chromosome %d, #Samples = %d' % (chrom, m.shape[0]), fontsize=18)
    P.xlabel('Sample Index C', fontsize=18)
    P.ylabel('Separation Measure M(C)', fontsize=18)
    i = np.argsort(m[:, 1])
    P.plot(m[i, 1], 'b-', label='M')
    P.ylim([-0.01, 1.01])
    # P.hold(True)
    # P.plot(m[i, 5], 'r-', label='Inbreeding')
    # P.legend(loc='lower right')
    P.show()
    P.savefig(os.environ['OBER'] + '/doc/poo/m-chr%d.png' % (chrom,))

def determine_poo(chrom, pedigree=None, genotyped_id_file=None, ibd_index_location=None, params=None, num_processes=1):
    '''Main call of determining parent-of-origin phases of quasi-founders, and flipping the samples
    that are found to be reversed.'''
    a = Aligner(chrom, pedigree=pedigree, genotyped_id_file=genotyped_id_file,
                ibd_index_location=ibd_index_location, params=params)
    return a.flip_measure_all_samples(threshold=params.poo_snp_measure_threshold, snp_step_size=params.poo_snp_step_size, debug=params.debug, num_processes=params.num_processes)

'''Dan''s kinship ratio function.'''
rat = lambda m: m[0] * m[3] / (m[1] * m[2])

####################################################################################
class Aligner(object):
    '''Determines the POO phase of quasi-founders for a single chromosome.'''

    '''Empty IBD clique (singleton).'''
    _EMPTY_CLIQUE = np.empty((0, 0))
    '''Kinship threshold for throwing outliers from the regression
    (if both father + mother kinships exceed this value, it's considered an outlier).''' 
    _K_THRESHOLD = 0.1  # 0.1 #0.05
    
    #---------------------------------------------
    # C-tor
    #---------------------------------------------
    def __init__(self, chrom, pedigree=None, genotyped_id_file=None, ibd_index_location=None, params=None):
        '''Load data and initialize for chromosome number chrom.'''
        self.params = params if params else im.param.PhaseParam()
        # TODO: replace this hard-coded value by pedigree location parameter passing
        self.ped = im.io_pedigree.read(pedigree, genotyped_id_file=genotyped_id_file) if pedigree else im.hutt_pedigree() 
        self.ibd = im.index.segment_index.SegmentIndex(ibd_index_location if ibd_index_location else os.environ['OBER_OUT'] + '/index_segments')
        self.num_snps = len(self.ibd._snp)
        # Quasi-founder = either parent of his/her is not genotyped 
        self.quasi_founders = np.where([all((y >= self.ped.num_genotyped) for y in self.ped.graph.predecessors(x))
                                        for x in xrange(self.ped.num_genotyped)])[0]
        self.quasi_founders_set = set(self.quasi_founders)
        self._init_chrom(chrom)

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def flip_measure_debug_sample(self, i, min_clique_size=3, threshold=0.25):
        '''Print debugging plots of flip measure at sample i.'''
        stats, k, snps = self.flip_measure_sample(self.ped, self.params, self.num_snps,
                                                  i, min_clique_size=min_clique_size, threshold=threshold)
        self.plot_separation_sample(stats, k, snps)
        return stats, k, snps

    def flip_measure_all_samples(self, min_clique_size=3, threshold=0.25, snp_step_size=1, samples=None,
                                 debug=False, num_processes=1):
        '''For each quasi-founder, generate a measure of flip (<1: need to flip; >1: no need to).
        Return an array of sample ids and corresponding flip measures.'''
        if num_processes == 1:  # Serial run
            return np.array([self.flip_measure_sample(self.ped, self.params, self.num_snps,
                                                      i, min_clique_size=min_clique_size, threshold=threshold,
                                                      snp_step_size=snp_step_size, debug=debug)[0]
                             for i in (samples if samples is not None else self.quasi_founders_set)])
        else:  # Parallel run
            raise ValueError('TBI')
#            manager = Manager()
#            lock = manager.Lock()
#            # Map phase
#            po = Pool(processes=num_processes)
#            m = po.imap(__determine_poo_sample, ((self.ped, self.params, sample, lock) for sample in self.quasi_founders))
#            m = np.array(m)
#            return m

    def flip_measure_sample(self, ped, params, num_snps, i, min_clique_size=3, threshold=0.25, snp_step_size=1, debug=False):
        '''Sample flip measure at all SNPs.'''
        snps = np.arange(0, num_snps, snp_step_size)
        # Dan recommends removing the proband's sibs that have a large kinship with both parents
        # and obstruct finding the right phase
        ignored_samples = ped.find_family_by_child(i, genotyped=False).children  # set([i])
        k = np.array(map(lambda x: (x[0],) + x[1],
                     (self.flip_measure(i, snp, min_clique_size=min_clique_size, ignored_samples=ignored_samples) for snp in snps)))
        ok = np.where((k[:, 1] >= min_clique_size) & (k[:, 2] >= min_clique_size))[0]
        k, snps = k[ok], snps[ok]
        sep = k[:, 0]
        num_snps = len(sep)
        n_plus, n_minus = len(np.where(sep > threshold)[0]), len(np.where(sep < -threshold)[0])
        m = max(n_plus, n_minus) / float(n_plus + n_minus + 1e-15)

        parents = [ped.sample_id[ped.graph.predecessors(i)[a]] for a in im.constants.ALLELES]
        inbreeding = params.kinship(parents[0], parents[1])

        stats = (i, m, n_plus, n_minus, num_snps, inbreeding, 1 if n_plus > n_minus else -1)
        if debug: print stats
        return stats, k, snps

    def ibd_clique_sample_ids(self, i, snp, a, ignored_samples=None):
        '''Return an array with the sample IDs of haplotypes in the IBD clique at sample i,
        SNP snp, allele a. Does not include i''s and i''s sibs IDs (or a custom non-None set of samples
        specified in ignored_samples).'''
        clique = self.ibd.find(self.chrom, snp, i, a)[:, 0]
        clique_without_i = list(set(clique) & self.quasi_founders_set - (ignored_samples if ignored_samples else self.ped.find_family_by_child(i, genotyped=False).children))
        return self.ped.sample_id[np.array(clique_without_i)] if clique_without_i else Aligner._EMPTY_CLIQUE
    
    def flip_measure(self, i, snp, min_clique_size=3, small_clique=10, ignored_samples=None):
        '''Suggested by Dan.
        Flip measure for sample i at a snp. This is the ratio between the two diagonals of the
        2x2 matrix of MEDIAN kinships (to exclude outliers with very large or small kinships with
        the parents). The ratio is scaled to a measure between [-1,1] (1=one phase, -1=opposite phase).'''
        c0, c1 = self.ibd_clique_sample_ids(i, snp, 0, ignored_samples=ignored_samples), self.ibd_clique_sample_ids(i, snp, 1, ignored_samples=ignored_samples)
        if min(c0.shape[0], c1.shape[0]) < min_clique_size:
            sep, k, c = None, None, None
        else:
            c = [c0, c1]
            parents = [self.ped.sample_id[self.ped.graph.predecessors(i)[a]] for a in im.constants.ALLELES]
            K = lambda k, l: np.array([self.params.kinship(parents[k], x) for x in c[l]])
            k = (K(0, 0), K(1, 0), K(0, 1), K(1, 1))
            sep = POSITIVE_TO_MINUS_1_TO_1(rat([np.median(x) for x in k]))
        return sep, (len(c0), len(c1)), (k, None, c)
    
    def flip_measure_regression(self, i, snp, min_clique_size=3, small_clique=10):
        '''Suggested by Oren.
        Flip measure for sample i at a snp. This is the separation angle between the regression slopes of
        the kinship coefficient of clique-father vs. the kinship coefficient of clique-mother for each
        of the two cliques. The angle is scaled to [-1,1].'''
        c0, c1 = self.ibd_clique_sample_ids(i, snp, 0), self.ibd_clique_sample_ids(i, snp, 1)
        if min(c0.shape[0], c1.shape[0]) < min_clique_size:
            sep, k, line, c = None, None, None, None
        else: 
            c = [c0, c1]
            parents = [self.ped.sample_id[self.ped.graph.predecessors(i)[a]] for a in im.constants.ALLELES]
            K = lambda k, l: np.array([self.params.kinship(parents[k], x) for x in c[l]])
            
            # Perform least-squares fit of both cliques
            # k_threshold = Aligner._K_THRESHOLD
            prob = [0.25, 0.75]  # Quantile region to include in fit            
            k00 = K(0, 0)
            k10 = K(1, 0)
            # good = (k00 <= k_threshold) | (k10 <= k_threshold)
            m0, m1 = mquantiles(k00, prob=prob), mquantiles(k10, prob=prob)
            good = (k00 >= m0[0]) & (k00 <= m0[1]) & (k10 >= m1[0]) & (k10 <= m1[1]) if len(c0) >= small_clique else np.arange(len(c0))
            a0, b0 = linefit(k00[good], k10[good], 1)
            # print len(c0), np.where(good)[0]
            
            k01 = K(0, 1)
            k11 = K(1, 1)
            # good = (k01 <= k_threshold) | (k11 <= k_threshold)
            m0, m1 = mquantiles(k01, prob=prob), mquantiles(k11, prob=prob)
            good = (k01 >= m0[0]) & (k01 <= m0[1]) & (k11 >= m1[0]) & (k11 <= m1[1]) if len(c1) >= small_clique else np.arange(len(c1))
            a1, b1 = linefit(k01[good], k11[good], 1)
            # print len(c1), np.where(good)[0]
            
            k = (k00, k10, k01, k11)
            line = (a0, b0, a1, b1)
            
            # Calculate the separation angle
            if np.abs(a0 - a1) < 1e-15:
                sep = 0
            else:
                e = (b1 - b0) / (a0 - a1)
                p, q = e, a0 * e + b0
                x0 = np.sign(a0 + b0)
                sep = max(-1, min(1, ANGLE_SCALING_FACTOR * (np.arctan2(a0 * x0 + b0 - q, x0 - p) - np.arctan2(a1 + b1 - q, 1 - p))))            
        return sep, (len(c0), len(c1)), (k, line, c)

    def plot_separation(self, i, snp, sep, data):
        '''Plot separation of IBD clique at sample i, SNP snp.''' 
        parents = [self.ped.sample_id[self.ped.graph.predecessors(i)[a]] for a in im.constants.ALLELES]
        print data
        print 'data[0]', data[0]
        print 'data[1]', data[1]
        (k00, k10, k01, k11), line = data[0:2]
    
        P.figure(1)
        P.clf()
        P.hold(True)
        # P.scatter([p], [q], color='k')
        P.scatter(k00, k10, color='r')
        P.scatter(k01, k11, color='b')
        if line is not None:
            a0, b0, a1, b1 = line        
            t = linspace(0, max(k00))
            P.plot(t, polyval([a0, b0], t), 'r')
            t = linspace(0, max(k01))
            P.plot(t, polyval([a1, b1], t), 'b')
        P.title('Sample %d, SNP %d, inbreeding=%.4f, separation=%.2f' % (i, snp, self.params.kinship(parents[0], parents[1]), sep))
        P.show()
        P.savefig(os.environ['OBER'] + '/doc/poo/sep-chr%d-%d-%d.png' % (self.chrom, i, snp))
        
    def plot_separation_sample(self, stats, k, snps):
        '''Plot separation of IBD clique at sample i, all SNPs along the chromosome.'''
        i, m = stats[0:2]
        P.figure(2)
        P.clf()
        P.hold(True)
        P.title('Scaled Separation Measure, Chromosome %d, Sample %d, $M$ = %.3f' % (self.chrom, i, m))
        P.xlabel('$s$ (SNP Index)')
        P.ylabel('$\phi(s)$')
        P.plot(snps, k[:, 0], 'b-')
        # P.plot(snps, np.minimum(k[:, 1], k[:, 2]), 'r-')
        P.show()
        P.savefig(os.environ['OBER'] + '/doc/poo/sep-angle-chr%d-%d.png' % (self.chrom, i))
        
    def debug_sample(self, i):
        stats, k, snps = self.flip_measure_sample(self.ped, self.params, self.num_snps, i)
        self.plot_separation_sample(stats, k, snps)
        return stats, k, snps
    
    def debug_sample_snp(self, i, snp):
        sep, _, data = self.flip_measure(i, snp)
        self.plot_separation(i, snp, sep, data)
        return snp, sep, data
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def _init_chrom(self, chrom):
        self.chrom = chrom
        self.total_snps = im.io.read_info_npz(os.environ['OBER_OUT'] + '/phasing/chr%d/hutt.phased.info.npz' % (chrom,)).num_snps
        self.ibd._load_chrom(chrom)
