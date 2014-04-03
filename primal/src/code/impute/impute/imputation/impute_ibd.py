#!/usr/bin/env python
'''
============================================================
Impute SNPs from 98 CGI WGS samples to all Hutt samples
using an IBD segment dictionary (a SmartSegmentSet).

Created on February 2, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, util, time, sys
from impute.data.constants import SNP, SAMPLE, MEGA_BASE_PAIR, ALLELES, INDETERMINATE
 
#---------------------------------------------
# Methods
#---------------------------------------------
def impute(h, ibd, t, snp=None, centrality_attenuation=0.9, debug=1, samples=None, ibd_sample_index=None):
    '''Main call that separately imputes each SNP in a list of SNPs.
    h = Haplotype object with phased haplotypes
    ibd = IBD dictionary - a SmartSegmentSet instance
    t = training genotypes. Must be compatible with problem''s pedigree
    samples = samples to impute (if None, all samples are imputed)
    centrality_attenuation = Attenuation factor to penalize non-central segments covering the SNP in question
    
    ibd_sample_index = optional dictionary of IBD-segment-ID-problem-sample-ID. If None, the mapping between the two IDs is assumed to be the identity.
    NOTE: NOT SUPPORTED YET - ONLY problem_to_ibd_id=None IS. TODO: ADD SUPPORT FOR GENERAL PROBLEMS. 
    '''
    if debug >= 1:
        print 'Imputing SNPs'
        print 'centrality_attenuation', centrality_attenuation
    snp = snp if snp is not None else np.arange(len(t.snp))
    # print 'snp', snp
    for snp_index in snp:
        t_start = time.time()
        chrom = t.snp['chrom'][snp_index]
        snp_bp = t.snp['base_pair'][snp_index]
        if debug >= 1:
            print '====== SNP %4d (%-22s): chr%-2d:%-9d x=%.2f ======' % \
        (snp_index, t.snp['name'][snp_index], chrom, snp_bp, t.snp['dist_cm'][snp_index])
        i = IbsImputer(h, ibd, t.training_data[snp_index, :, :],
                       t.sample_index, snp_bp, centrality_attenuation, debug=(debug >= 2))
        result = i.result
        result = i.impute(samples)
        t.imputed_data[snp_index, :, :] = result
        if debug >= 2:
            np.set_printoptions(threshold=np.nan)
            print list(enumerate(result.tolist()))
        if debug >= 1:
            print 'Time', time.time() - t_start
            print 'Overall allele call rate %.2f%%' % ((100.0 * len(result.nonzero()[0])) / result.size,)
            print 'Imputed', len(result.nonzero()[0]), 'out of', result.size 
#        data = h.data[snp_index, :, :]
#        c1 = len(np.where(data == 1)[0])
#        c2 = len(np.where(data == 2)[0])
#        print 'Allele frequencies %.2f, %.2f' % ((1.0 * c1) / data.size, (1.0 * c2) / data.size)
    return i

####################################################################################
class IbsImputer(object):
    '''Imputes a single base-pair location using IBD segments.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, h, ibd, g, training_sample_index, snp_bp, centrality_attenuation, debug,
                 num_passes_training=2, num_passes_target=1):
        '''Initialize an imputer from phasing data h and training genotype data g. Imputing
        location h[k,:,:], which is assumed to correspond to g for the training data.
        
        num_passes_target>1 has an effect only if all 1415x1415 segments are used. In principle,
        both passes should be set to 1 if there is perfect IBD transitivity. In practice, we have
        near-transitivity, so 2 passes instead of 1 can be beneficial.'''
        # One underscore = protected member
        self.h = h
        self.g = g
        self.ibd = ibd
        self.training_sample_index = training_sample_index
        self.snp_bp = snp_bp
        self.snp = np.arange(h.data.shape[0])
        self.num_passes_training = num_passes_training
        self.num_passes_target = num_passes_target

        # Store the imputation result in an array that looks like one SNP in h
        self.result = np.zeros_like(h.data[0, :, :])
        
        # Maps sample ID to training set index
        self.training_index = dict(zip(self.training_sample_index, xrange(self.g.shape[0])))
        self.bp = self.h.snp['base_pair']
        # self.stats = []
        self.__mu = centrality_attenuation
        self.debug = debug
         
    #---------------------------------------------
    # Methods
    #---------------------------------------------        
    def impute(self, samples=None):
        '''Infer imputed genotypes at all samples of h from the samples of g, at the location
        between k-1 and k.'''
        
        # Phase all hom training samples         
        self.__phase_hom()

        # Phase as much as possible non-hom training samples
        # Bootstrap: impute R -> T; pass imputed samples from T to R; repeat
        R = azip(*util.flattened_meshgrid(self.hom, ALLELES))
        T = set(self.non_hom)
        R, _ = self._impute_bootstrap(R, T, self.num_passes_training, True, 'non-hom typed')
        
        # Impute target  samples from all available training haps
        T = set(range(self.h.num_samples)) - set(self.training_sample_index)
        if samples is not None:
            T &= set(samples)
        R, T = self._impute_bootstrap(R, T, self.num_passes_training, False, 'non-typed')

        return self.result
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def _impute_bootstrap(self, R, T, num_passes, genotype_available, title):
        '''Impute the target sample set T using the R-haps. Perform num_passes bootstrap passes.
        Return the updated R and T (every imputed hap is transferred from T to R).'''
        for i in xrange(num_passes):
            if self.debug:
                print '#' * 80
                print 'Imputing %s pass %d/%d: |R| %d |T| %d' % \
                (title, i + 1, self.num_passes_training, len(R), len(T))
                print '#' * 80
                print 'R', R
                print 'T', T
            # Gauss-Seidel-type pass: add any imputed sample to R within the pass
            imputed = []
            for sample in T:
                a = self._impute_sample(sample, R, genotype_available)
                if a:
                    imputed.append(sample)
                    R += a
            T -= set(imputed)
        return R, T

    def __phase_hom(self):
        '''Phase all homozygous training samples; place in corresponding locations in h = imputed_data.
        Set the training sample ID subsets hom, non_hom.'''
        is_hom = im.gt.is_homozygous(self.g)[:]
        hom_index = np.where(is_hom)
        self.hom = self.training_sample_index[hom_index]
        self.non_hom = self.training_sample_index[np.where(~is_hom)]
        self.result[self.hom, :] = self.g[hom_index]

    def _impute_sample(self, sample, R, genotype_available):
        '''Phase and impute non-homozygous training sample from a phased training hap set R,
        if genotype_available=True; otherwise, impute a non-training sample from R.'''
        if self.debug:
            print '*' * 15
            print 'sample', sample
            print '*' * 15
        imputed_alleles = []
        hap_nbhr = [INDETERMINATE, INDETERMINATE]
        hap_L = [INDETERMINATE, INDETERMINATE]
        
        # For each allele, find the longest IBD segment and corresponding neighbor index in R
        for a in ALLELES:
            nbhr, L = self._ibd_sharing_sample((sample, a), R)
            hap_nbhr[a], hap_L[a] = nbhr, L
            if self.debug:
                if not L:
                    sys.stdout.write('No segments found\n')
                else:
                    # self.stats.append(((sample, a), nbhr, R[SNP][nbhr], R[SAMPLE][nbhr], L, nbhr != INDETERMINATE))
                    nbhr_hap = self.result[nbhr[SNP], nbhr[SAMPLE]]
                    sys.stdout.write('(%-4d, %1d) [%-4d, %1d) L %7.4f Mbp %s' % \
                        (sample, a, nbhr[SNP], nbhr[SAMPLE], L / MEGA_BASE_PAIR, 'YES' if nbhr != INDETERMINATE else '-'))
                    if nbhr != INDETERMINATE:
                        sys.stdout.write('\tallele: %d' % nbhr_hap)
                    sys.stdout.write('\n')
        
        # Impute; if phasing a genotyped sample, use the longer IBD segment of the two alleles
        for a in ((0, 1) if hap_L[0] > hap_L[1] else (1, 0)):
            nbhr, L = hap_nbhr[a], hap_L[a]  
            if nbhr != INDETERMINATE:
                # Found IBD-sharing nbhr hap
                # Phase corresponding entry
                self.result[sample, a] = self.result[nbhr[SNP], nbhr[SAMPLE]]
                imputed_alleles.append((sample, a))
                if self.debug:
                    sys.stdout.write('Updated hap entry (%-4d, %1d) = %d%d based on L=%7.4f Mbp, nbhr = (%-4d, %1d)\n' % \
                                     (sample, a, self.result[sample, 0], self.result[sample, 1], L / MEGA_BASE_PAIR, nbhr[SNP], nbhr[SAMPLE]))
                if genotype_available:
                    # Impute other entry if we have a second genotype
                    im.gt.complete_haplotype_single(self.result[sample, :],
                                                    self.g[self.training_index[sample], :])
                    imputed_alleles.append((sample, 1 - a))
                    if self.debug:
                        sys.stdout.write('Updated hap entry (%-4d, %1d) = %d%d based on L=%7.4f Mbp, nbhr = (%-4d, %1d)\n' % \
                                         (sample, a, self.result[sample, 0], self.result[sample, 1], L / MEGA_BASE_PAIR, nbhr[SNP], nbhr[SAMPLE]))
                    break
        # Returns list of imputed haplotype indices
        return imputed_alleles

    def _segment_len(self, segments):
        '''Return segment endpoint and length arrays.'''
        a = np.array(map(lambda x: x.start, segments))
        b = np.array(map(lambda x: x.stop, segments))
        return a, b, b - a

    def _segment_attn(self, a, b, L):
        '''Attenuate segment size by a factor of mu if k is at the edge of the interval
        (no penalty in the middle), so that central segments are preferred.''' 
        return L * (self.__mu + (1 - self.__mu) * np.abs(self.snp_bp - 0.5 * (a + b)) / (0.5 * (b - a)))

    def _print_segment_scores(self, a, b, L, segment_attn_len, segments):
        '''Print segment scores.'''
        print a, b, L, segment_attn_len
        stats = np.concatenate((
                                a[np.newaxis].transpose(),
                                b[np.newaxis].transpose(),
                                L[np.newaxis].transpose(),
                                segment_attn_len[np.newaxis].transpose(),
                                np.array(map(lambda x: x.samples[1][0], segments))[np.newaxis].transpose(),
                                np.array(map(lambda x: x.samples[1][1], segments))[np.newaxis].transpose()
                                ), axis=1)
        stats = stats[np.lexsort((stats[:, 0], stats[:, 1], -stats[:, 3]))]
        # Print longest segments' information
        for x in stats[0:10]:
            print '\t(%-4d,%1d) [%-4d, %1d) %5.2f %-5.2f (%d)' % \
            (x[4], x[5], x[0], x[1], x[2] / MEGA_BASE_PAIR, x[3] / MEGA_BASE_PAIR, self.result[int(x[4]), int(x[5])])

    def _ibd_sharing_sample(self, hap, R):
        '''Return the index of an IBD-sharing sample with the haplotype hap=(sample,allele) among
        a list of phased training samples R = (array-of-samples, array-of-alleles),
        or MISSING, if none found or the IBD segment [bp] < min_ibs_len_mbp.'''
        # Find all segments between hap and R that intersect the SNP position
#        intersecting = self.ibd.find(self.snp_bp, self.snp_bp)
#        segments = self._relevant_intersecting(intersecting, hap, R)
        segments = self.ibd.find(self.snp_bp, self.snp_bp, hap, R) 

        # Find the best IBD segment and return its index if it's sufficiently long
        a, b, L = self._segment_len(segments)    
        if not L.size:
            return INDETERMINATE, None
        
        # Attenuate segment size by a factor of mu if k is at the edge of the interval
        # (no penalty in the middle), so that central segments are preferred 
        segment_attn_len = self._segment_attn(a, b, L)
        nbhr = np.argmax(segment_attn_len)
        if self.debug:
            # print 'Intersecting segments', len(intersecting)
            print 'Relevant Segments', len(segments)
            self._print_segment_scores(a, b, L, segment_attn_len, segments)
        return segments[nbhr].samples[1], L[nbhr]

def azip(*args):
    '''Zip a list of arrays.'''
    return zip(*(x.tolist() for x in args))
