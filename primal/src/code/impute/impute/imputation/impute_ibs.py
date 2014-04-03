#!/usr/bin/env python
'''
============================================================
Impute SNPs from 98 CGI WGS samples to all Hutt samples
using long identity-by-state segments. A workaround
suggested by Dan to address the immediate need for Carole's
grant proposal on rare variants due 15-DEC-2012.

Created on December 3, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, util, time, sys
from impute.data.constants import SNP, SAMPLE, MEGA_BASE_PAIR, ALLELES, INDETERMINATE
from scipy import ndimage
from utilities.math.index_util import first_occurrence_index_byte
 
#---------------------------------------------
# Methods
#---------------------------------------------
def impute(reader, t, snp=None, min_ibs_len_mbp=1.5, max_ibs_len_snp=400, error_filter_length=5,
           centrality_attenuation=0.9, debug=False):
    '''Main call that separately imputes each SNP in a list of SNPs.
    min_ibs_len_mbp = Minimum IBS segment required for imputation [Mbp]
    max_ibs_len_snp = Largest IBS segment of interest [SNPs] = size of window to search for IBS segments around imputed snp
    error_filter_length = Median filter size to ignore isolated errors and non-phased haplotypes
    centrality_attenuation = Attenuation factor to penalize non-central segments covering the SNP in question 
    '''
    print 'Imputing SNPs'
    print 'min_ibs_len_mbp', min_ibs_len_mbp, 'max_ibs_len_snp', max_ibs_len_snp, 'error_filter_length', error_filter_length
    snp = snp if snp is not None else np.arange(len(t.snp))
    # print 'snp', snp
    for snp_index in snp:
        t_start = time.time()
        chrom = t.snp['chrom'][snp_index]
        bp = t.snp['base_pair'][snp_index]
        print '====== SNP %4d (%-22s): chr%-2d:%-9d ======' % (snp_index, t.snp['name'][snp_index], chrom, bp)
        
        # Get a window around the SNP wide enough to accommodate the largest IBS segment of interest
        h, k = reader.get_window(chrom, bp, max_ibs_len_snp)
        print 'num_snps', h.num_snps, 'relative k', k
        
        i = IbsImputer(h, t.training_data[snp_index, :, :], t.sample_index, k,
                       min_ibs_len_mbp, error_filter_length, centrality_attenuation, debug=debug)
        result = i.impute()
        t.imputed_data[snp_index, :, :] = result
        print 'Time', time.time() - t_start
        print 'Overall allele call rate %.2f%%' % ((100.0 * len(result.nonzero()[0])) / result.size,)
    return i

####################################################################################
class IbsImputer(object):
    '''Imputes a single base-pair location using IBD segments.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, h, g, training_sample_index, k, min_ibs_len_mbp, error_filter_length,
                 centrality_attenuation, debug):
        '''Initialize an imputer from phasing data h and training genotype data g. Imputing
        location h[k,:,:], which is assumed to correspond to g for the training data.'''
        # One underscore = protected member
        self.h = h
        self.g = g
        self.training_sample_index = training_sample_index
        self.k = k
        self.min_ibs_len_mbp = MEGA_BASE_PAIR * min_ibs_len_mbp
        self.error_filter_length = error_filter_length
        self.snp = np.arange(h.data.shape[0])

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
    def impute(self):
        '''Infer imputed genotypes at all samples of h from the samples of g, at the location
        between k-1 and k.'''
        
        # Phase all hom training samples         
        # print '#' * 80
        
        # print '#' * 80
        self.__phase_hom()
        print 'Imputing genotyped hom', '|hom|', len(self.hom), '|non_hom|', len(self.non_hom)
        
        # Phase as much as possible non-hom training samples
        R = util.flattened_meshgrid(self.hom, ALLELES)
        # set_printoptions(threshold=np.nan)
        # R = [np.array([1053]), np.array([0])]
        if self.debug:
            print '#' * 80
        print 'Imputing genotyped non-hom', '|R|', len(R[0]), '|T|', len(self.non_hom)
        if self.debug:
            print '#' * 80
            print 'R', R
            print 'T', self.non_hom
        for sample in self.non_hom:
            self._impute_sample(sample, R, True)
        
        # Impute non-training samples from all available training haps
        # R = [np.array([1053, 0]), np.array([0])]
        R = self.result.nonzero()
        T = set(range(self.h.num_samples)) - set(self.training_sample_index)
        if self.debug:
            print '#' * 80
        print 'Imputing non-genotyped', '|R|', len(R[0]), '|T|', len(T)
        if self.debug:
            print '#' * 80
            print 'R', R
            print 'T', self.non_hom
        for sample in T:
            self._impute_sample(sample, R, False)
        return self.result
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------        
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
        hap_nbhr = [INDETERMINATE, INDETERMINATE]
        hap_L = [INDETERMINATE, INDETERMINATE]
        
        # For each allele, find the longest IBS segment and corresponding neighbor index in R
        for a in ALLELES:
            nbhr, L = self._ibs_sharing_sample((sample, a), R)
            nbhr = nbhr if L >= self.min_ibs_len_mbp else INDETERMINATE
            hap_nbhr[a], hap_L[a] = nbhr, L
            if self.debug:
                #self.stats.append(((sample, a), nbhr, R[SNP][nbhr], R[SAMPLE][nbhr], L, nbhr != INDETERMINATE))
                nbhr_hap = self.result[R[SNP][nbhr], R[SAMPLE][nbhr]]
                sys.stdout.write('(%-4d, %1d) nbhr %3d (%-4d, %1d) L %7.4f Mbp %s' % \
                    (sample, a, nbhr, R[SNP][nbhr], R[SAMPLE][nbhr], L / MEGA_BASE_PAIR, 'YES' if nbhr != INDETERMINATE else '-'))
                if nbhr != INDETERMINATE:
                    sys.stdout.write('\tallele: %d' % nbhr_hap)
                sys.stdout.write('\n')
        
        # Impute; if phasing a genotyped sample, use the longer IBS segment of the two alleles
        for a in ((0, 1) if hap_L[0] > hap_L[1] else (1, 0)):
            nbhr, L = hap_nbhr[a], hap_L[a]  
            if nbhr != INDETERMINATE:
                # Found IBS-sharing nbhr hap
                # Phase corresponding entry
                self.result[sample, a] = self.result[R[SNP][nbhr], R[SAMPLE][nbhr]]
                if self.debug:
                    sys.stdout.write('Updated hap entry (%-4d, %1d) = %d%d based on L=%7.4f Mbp\n' % \
                                     (sample, a, self.result[sample, 0], self.result[sample, 1], L / MEGA_BASE_PAIR))
                if genotype_available:
                    # Impute other entry if we have a second genotype
                    im.gt.complete_haplotype_single(self.result[sample, :],
                                                    self.g[self.training_index[sample], :])
                    if self.debug:
                        sys.stdout.write('Updated hap entry (%-4d, %1d) = %d%d based on L=%7.4f Mbp\n' % \
                                         (sample, a, self.result[sample, 0], self.result[sample, 1], L / MEGA_BASE_PAIR))
                    break

    def _ibs_sharing_sample(self, hap, R):
        '''Return the index of an IBS-sharing sample with the haplotype hap=(sample,allele) among
        a list of phased training samples R = (array-of-samples, array-of-alleles),
        or MISSING, if none found or the IBS segment [bp] < min_ibs_len_mbp.'''
        
        # Calculate difference between hap and each element in R
        h_sample = self.h.data[:, hap[SNP], hap[SAMPLE]]  # Haplotype of target sample 
        h_r = self.h.data[:, R[SNP], R[SAMPLE]]  # Haplotypes of phased training samples
        d = np.transpose(np.tile(h_sample, (h_r.shape[1], 1)).transpose() != h_r).astype('int8')
        # Ignore unphased SNPs
        # filled = np.where(d != INDETERMINATE) # May be used in the future to determine size of segments in snp s
        unphased = np.where(d == INDETERMINATE)
        d[unphased] = 0
        
        # Filter to ignore genotype errors and isolated missing phasing
        # (note: abs value function only needed to map MISSING, errors to 1 if they are not filtered above)
        d = ndimage.median_filter(d, 'median', footprint=np.ones((1, self.error_filter_length)))
        # d = ndimage.median_filter(abs(d), 'median', footprint=np.ones((1, self.error_filter_length)))
        # for i in xrange(d.shape[0]):
            # d[i, :] = filter_diff(d[i, :], 'fast_max', self.error_filter_length)
        
        # Find the longest IBS segment and return its index if it's sufficiently long    
        segment = np.array([self._ibs_segment_len(dr, self.k) for dr in d])
        
        # Attenuate segment size by a factor of mu if k is at the edge of the interval
        # (no penalty in the middle), so that central segments are preferred 
        (a, b, L) = (segment[:, 0], segment[:, 1], segment[:, 2])
        all_segments = L * (self.__mu + (1 - self.__mu) * np.abs(self.k - 0.5 * (a + b)) / (0.5 * (b - a)))
        nbhr = np.argmax(all_segments)
        if self.debug:
            stats = np.concatenate((all_segments[np.newaxis].transpose(), R[SNP][np.newaxis].transpose(), R[SAMPLE][np.newaxis].transpose()), axis=1)
            stats = stats[np.lexsort((stats[:, 2], stats[:, 1], -stats[:, 0]))]
            # Print longest segments' information
            for x in stats[0:10]:
                print '\t(%-4d, %1d) %5.2f' % (x[1], x[2], x[0] / MEGA_BASE_PAIR)
        return nbhr, L[nbhr]
    
    def _ibs_segment_len(self, d, k):
        '''Find the largest segment of consecutive zeros containing k. Returning the
        length in bp. snp is the snp index array.'''
        left = first_occurrence_index_byte(d, 1, k - 1, -1)
        left = left if left >= 0 else 0
        right = first_occurrence_index_byte(d, 1, k, 1)
        right = right if right >= 0 else len(d) - 1
        L = self.bp[right] - self.bp[left]
        #if self.debug:
        #    print '\t', 'k', k, 'd around k', d[k - 10:k + 10]
        #    print '\t', left, right, self.bp[left], self.bp[right], L
        return left, right, L
