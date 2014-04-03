'''
============================================================
Generalized Mendelian error check of the WGS CGI data
Quality Control (QC). Uses cliques to assign an error rate
(or concordance) to each variant found in the 98 sequenced
Hutterites.

Similar to imputation (cf. impute_ibd_index.py), but only
uses the data on the 98. All alleles are attempted to be
assigned to a clique. The number of alleles who do not agree
with a clique's consensus are counted as errors.

Created on November 6, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, time, itertools, os
from impute.data.constants import MISSING, PATERNAL, MATERNAL, ALLELES, UNPHASED
from Queue import Queue
from impute.tools.genotype_tools import complete_haplotype_partial
from impute.tools import recode
from collections import Counter
from impute.imputation.reader import problem_to_imputation_set
from impute.data.problem import Problem
from db_gene.cgi.ids import DEFAULT_ID_INDEX
from impute.phasing.examples import OUT_PHASING
 
# Indices into the dimension of a haplotype array of a single SNP
SAMPLE, ALLELE = range(2)

#---------------------------------------------
# Methods
#---------------------------------------------
def qc_problem(problem, index_file=DEFAULT_ID_INDEX, samples=None, snps=None, debug=0,
               input_dir=None, segment_location=None, ibd_sample_index=None, debug_sample= -1):
    '''Impute the SNPs that have also been genotyped for all samples, and compare with the
    imputation result. Returns the phasing Problem and the imputed results.
    
    problem = npz file location or Problem instance.
    segment_location = IBD segment repository (segment file or segment index) location.
    ibd_sample_index = optional dictionary of IBD-segment-ID-problem-sample-ID. If None, the mapping between the two IDs is assumed to be the identity.
    
    algorithm options:
    'index' - Imputation V2 using IBD cliques at each SNP (default).
    'interval_tree' - Imputation V1 using interval tree queries on a SmartSegmentSet.'''
    # Read input data
    p = problem if isinstance(problem, Problem) else im.io.read_npz(problem)
    
    # Select SNPs
    if snps is not None: p = p.sub_problem_of_snps(snps)
    t = problem_to_imputation_set(p, index_file=index_file)
    print p
    print t.genotype

    # Load IBD index    
    ibd = im.index.segment_index.SegmentIndex(segment_location if segment_location else '%s/index_segments' % (os.environ['OBER_OUT'],))
    # Run imputation
    qc(ibd, t, samples=samples, debug=debug, genotype=p.g, ibd_sample_index=ibd_sample_index, debug_sample=debug_sample)
    return p, t

def qc_chrom(chrom=None, snp_step_size=None, snps=None, debug=0, input_dir=None, segment_location=None):
    '''Impute a phased chromosome. Uses default file locations. Imputed SNPs that are uniformly
    spaced (in snp index) with step size snp_step_size.
    Returns Problem, ImputationSet, snp range objects.'''
    input_dir = input_dir if input_dir else '%s/chr%d' % (OUT_PHASING, chrom)
    problem = im.io.read_npz('%s/hutt.phased.npz' % (input_dir,))
    snps = snps if snps is not None else (np.arange(0, problem.num_snps, snp_step_size) if snp_step_size else None)
    return qc_problem(problem, input_dir=input_dir, segment_location=segment_location, snps=snps, debug=debug)

def qc(ibd, t, snp=None, majority_threshold=0.66, debug=1, samples=None, genotype=None,
       debug_sample= -1, ibd_sample_index=None):
    '''Main call that separately imputes each SNP in a list of SNPs.
    ibd = IBD dictionary - a SegmentIndex instance
    t = ImputationSet instance. Contains training genotypes and imputed genotypes.
    majority_threshold = Threshold for resolving IBD conflicts via the majority vote 
    samples = samples to impute (if None, all samples are imputed)
    genotype = true genotype at all samples (non-None for validation studies).'''
    # If problem sample indices require translation, augment the t.imputed_data, t.imputed_hap_type arrays
    # to the #genotyped samples in the IBD segment index; then impute all of them; in the end, restrict the
    # result to the problem's sample index subset.
    imputed_data, imputed_hap_type = t.imputed_data, t.imputed_hap_type
    sample_translation = ibd_sample_index is not None
    if sample_translation:
        imputed_data = np.zeros((t.num_snps, ibd.num_samples, 2), dtype=np.byte)
        imputed_data[:, ibd_sample_index, :] = t.imputed_data
        imputed_hap_type = np.zeros((t.num_snps, ibd.num_samples), dtype=np.byte)
        imputed_hap_type[:, ibd_sample_index] = t.imputed_hap_type
        
    snp = snp if snp is not None else np.arange(len(t.snp))
    if debug >= 1:
        print 'Imputing SNPs'
        print 'majority_threshold', majority_threshold
        print 'snp', snp
    for snp_index in snp:
        chrom = t.snp['chrom'][snp_index]
        snp_bp = t.snp['base_pair'][snp_index]
        if debug >= 1:
            print '====== SNP %4d (%-22s): chr%-2d:%-9d x=%.2f ======' % \
        (snp_index, t.snp['name'][snp_index], chrom, snp_bp, t.snp['dist_cm'][snp_index])
        t_start = time.time()
        _IbdQc(imputed_data[snp_index], imputed_hap_type[snp_index],
               ibd, t.training_data[snp_index], t.sample_index,
               chrom, snp_bp, debug=(debug >= 2), majority_threshold=majority_threshold,
               debug_sample=debug_sample).impute(samples)           
        t_impute = time.time() - t_start
        if debug >= 1:
            result = imputed_data[snp_index]
            r = recode.recode_single_genotype(result)
            if sample_translation: r = r[ibd_sample_index]
            num_fully_called = len(np.where(r > 0)[0])
            num_alleles_called = len(result.nonzero()[0])
            result_training = result[t.sample_index]
            num_phased = len(result_training.nonzero()[0])
            # Per Rebecca's request: print the list of samples IDs that were not fully called
#             snp_name = t.snp['name'][snp_index]
#             for sample in im.examples.sample_index_to_id()[np.where(r <= 0)[0]]:
#                 print '%s,%s' % (snp_name, sample)
            if genotype is not None:
                rg = recode.recode_single_genotype(genotype[snp_index])
                if debug >= 2:
                    eq = ((r <= 0) | (rg <= 0) | (r == rg)).astype(int)  
                    np.set_printoptions(threshold=np.nan)
                    group_index = ibd._group_index[0 if ibd.test_index else ibd.nearest_left_snp(chrom, snp_bp) - ibd._start]
                    if sample_translation: group_index, result = group_index[ibd_sample_index], result[ibd_sample_index]
                    all_haps = np.concatenate((np.arange(t.num_samples)[np.newaxis].transpose(),
                                          genotype[snp_index], result.tolist(),
                                          group_index  # ,eq[np.newaxis].transpose(
                                          ), axis=1) 
                    print 'All (true vs. imputed):'
                    print all_haps
                    print 'Discordant (true vs. imputed):'
                    discordant = all_haps[np.where(eq == 0)[0]]
                    print discordant
                    bad_paternal = np.where(discordant[:, 1] != discordant[:, 3])[0]
                    bad_maternal = np.where(discordant[:, 2] != discordant[:, 4])[0]
                    print 'bad_paternal', bad_paternal
                    print 'bad_maternal', bad_maternal
    
                    # dis = compressed form of haplotype report obtained from the discordant array 
                    dis = np.zeros((discordant.shape[0], 3), dtype=np.uint)
                    if bad_paternal.size:
                        dis[bad_paternal, :] = discordant[bad_paternal][:, np.array([0, 1, 5])]
                    if bad_maternal.size:
                        dis[bad_maternal, :] = discordant[bad_maternal][:, np.array([0, 2, 6])]
                    dis = np.concatenate((dis, np.zeros((dis.shape[0], 1), dtype=np.uint)), axis=1)
                    if bad_paternal.size:
                        dis[bad_paternal, 3] = PATERNAL
                    dis[bad_maternal, 3] = MATERNAL
                    dis = dis[:, np.array([0, 3, 1, 2])]
                    bad_groups = np.array(Counter(dis[:, 3]).items(), dtype=np.uint)
                    if bad_groups.size:
                        bad_groups = bad_groups[np.lexsort((bad_groups[:, 0], -bad_groups[:, 1]))]
                    print 'bad_group id, #discordances\n', bad_groups

            G, H = t.training_data[snp_index], result[t.sample_index]
            changed = np.where((H[:, 0] != MISSING) & (H[:, 1] != MISSING) & (G[:, 0] != MISSING) & (G[:, 1] != MISSING) & 
                               (H[:, 0] + H[:, 1] != G[:, 0] + G[:, 1]))[0]

            print 'Time               %6.2f s' % (t_impute,)
            print 'Call rate allele   %6.2f%% (%d/%d)' % \
            ((100.0 * num_alleles_called) / result.size, num_alleles_called, result.size)
            print 'Call rate genotype %6.2f%% (%d/%d)' % \
            ((100.0 * num_fully_called) / result.shape[0], num_fully_called, result.shape[0])
            print 'Phased training    %6.2f%% (%d/%d)' % ((100.*num_phased) / result_training.size, num_phased, result_training.size)                
            if genotype is not None:
                c, con, dis = recode.concordance_recoded(r, rg)
                print 'Concordance        %6.2f%% (%d/%d)' % (100.*c, con, con + dis)
            if changed.size:
                print 'Changed training   %6.2f%% (%d/%d) %s' % ((100.*len(changed)) / len(t.sample_index), len(changed), len(t.sample_index), repr(t.sample_index[changed])[6:-1])
    
    # Restrict imputed results to problem's sample index subset
    if sample_translation:
        t.imputed_data = imputed_data[:, ibd_sample_index, :]
        t.imputed_hap_type = imputed_hap_type[:, ibd_sample_index]
        
####################################################################################
class _IbdQc(object):
    '''Calculates a QC measure for a single variant using IBD cliques.'''
    
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    __EMPTY_ARRAY = np.array([])
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, h, hap_type, ibd, g, training_sample_index, chrom, snp_bp, debug=False, majority_threshold=0.66,
                 debug_sample= -1, max_iter=1000):
        '''Initialize an imputer that changes the result h in-place using the IBD index ibd
        and training genotype data g at snp position snp_bp. majority vote = threshold for majority vote. When
        |# haps with majority allele| >= majority_threshold*(All haps) in a clique, the vote is accepted.''' 

        # Input fields
        self.h, self.hap_type, self.g, self.ibd, self.training_sample_index, self.max_iter, self.debug, \
        self.debug_sample, self.chrom = h, hap_type, g, ibd, training_sample_index, max_iter, debug, \
        debug_sample, chrom

        # Maps sample ID to training set index
        self.training_index = dict(zip(self.training_sample_index, xrange(self.g.shape[0])))
        self.ratio_threshold = majority_threshold / (1 - majority_threshold)

        # Find the appropriate IBD index SNP for the target base-pair position
        self.snp = ibd.nearest_left_snp(chrom, snp_bp)
        if self.debug:
            ibd.find(self.chrom, self.snp, self.training_sample_index[0], PATERNAL)
            print 'IBD index file %s/chr%d/region-%d.npz, nearest SNP %d' % (ibd._index_dir, ibd._chrom, ibd._start, self.snp)
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------                 
    def impute(self, samples=None):
        '''Infer imputed genotypes at all samples of h from the samples of g.'''
        # Aliases
        r = self.ratio_threshold
        if self.debug_sample >= 0:
            j = self.debug_sample
            if j in self.training_index: print 'Initial g[%d] = %s, h[%d] = %s' % (self.training_index[j], repr(self.g[self.training_index[j]]), j, repr(self.h[j]))
            else: print 'Initial h[%d] = %s' % (j, repr(self.h[j]))
        
        # Create a queue of haplotypes that can be used to impute others. Initially, it
        # is the training set homozygotes; every time we phase a training set het, its other
        # allele is appended to the queue.
        #
        # TODO: possibly replace by a priority queue where alleles are ordered by their clique sizes?
        # (we have extra confidence in those alleles; not sure it matters though) 
        q = Queue()

        # Initial condition: phase all hom training samples         
        hom = self.__phase_hom()
        for hap in itertools.product(hom, ALLELES):
            if self.debug: print 'Adding hom haplotype to queue', hap
            q.put(hap)
        num_hom_haps = q.qsize()
        if self.debug:
            print 'Items on queue : %d' % (q.qsize(),)
            print 'filled haps    : %.2f%%' % ((100.*len(self.h.nonzero()[0])) / self.h.size,)
            HH = recode.recode_single_genotype(self.h)
            print 'filled samples : %.2f%%' % ((100.*len(np.where(HH > 0)[0])) / HH.size,)
            print 'phased training: %.2f%%' % ((100.*len(self.h[self.training_sample_index].nonzero()[0])) / self.h[self.training_sample_index].size,)
        count = 0
        while not q.empty():
            # Find group = the IBD clique of an imputed haplotype (all haps that are IBD with it at self.snp)
            hap = q.get()
            count += 1
            if count > self.max_iter: raise ValueError('Possible imputation bug - exceeded maximum number of iterations')
            if self.debug:
                print '*' * 55
                print 'Iteration %d, imputing from hap %s' % (count, hap)
                print '*' * 55
                if self.debug_sample >= 0:
                    j = self.debug_sample
                    print 'h[%d] = %s' % (j, repr(self.h[j]))
                    if self.h[j, 0] == 1:
                        pass
            group = self.ibd.find(self.chrom, self.snp, hap[SAMPLE], hap[ALLELE])
            if group.size:
                # if self.debug: print 'group', group
                s, a = group[:, SAMPLE], group[:, ALLELE]
                H = self.h[s, a]
#                    print 'H', H
                # Find haplotypes that have been imputed as allele 1 and those imputed as allele 2
                R1, R2 = group[np.where(H == 1)], group[np.where(H == 2)]
                if self.debug:
                    print 'IBD group %d (%d haps):\n%s' % (self.ibd.group_index(self.chrom, self.snp, hap[SAMPLE], hap[ALLELE]), len(group), repr(np.concatenate((group, H[np.newaxis].transpose()), axis=1)))
                    print 'R1 = %s' % (repr(list(map(tuple, R1))),)
                    print 'R2 = %s' % (repr(list(map(tuple, R2))),)
#                    print 'R1:\n%s' % (repr(np.concatenate((R1, self.h[R1[:, 0], R1[:, 1]][np.newaxis].transpose()), axis=1)))
#                    print 'R2:\n%s' % (repr(np.concatenate((R2, self.h[R2[:, 0], R2[:, 1]][np.newaxis].transpose()), axis=1)))
                # Majority vote: if there are enough haps with one allele, override the rest.
                # Otherwise, an unresolved conflict ==> zero everyone out.
                l1, l2 = len(R1), len(R2)
                consensus = 1 if l1 >= r * l2 else (2 if l2 >= r * l1 else MISSING)
                self.h[s, a] = consensus
                if consensus == 0:
                    # If no consensus is reached, keep the already-imputed values in place, otherwise
                    # we can run into an infinite loop by imputing and erasing h-entries.  
                    self.h[R1[:, 0], R1[:, 1]] = 1
                    self.h[R2[:, 0], R2[:, 1]] = 2
                H = self.h[s]
                if self.debug:
                    print 'l1 %d l2 %d consensus %d' % (l1, l2, consensus)
                    print 'Items on queue : %d' % (q.qsize(),)
                    print 'filled haps    : %.2f%%' % ((100.*len(self.h.nonzero()[0])) / self.h.size,)
                    HH = recode.recode_single_genotype(self.h)
                    print 'filled samples : %.2f%%' % ((100.*len(np.where(HH > 0)[0])) / HH.size,)
                    print 'phased training: %.2f%%' % ((100.*len(self.h[self.training_sample_index].nonzero()[0])) / self.h[self.training_sample_index].size,)
                
                # Phase training hets (this includes BOTH partially-called = potential hets an
                # fully-called hets) with one imputed allele
                i = np.array([self.training_index.has_key(x) for x in s])
                si = s[i]
                G = self.g[map(self.training_index.get, si), :]
                unphased_hets = np.where(((H[i, PATERNAL] != MISSING) ^ (H[i, MATERNAL] != MISSING)) 
                                         & ((G[:, PATERNAL] != MISSING) | (G[:, MATERNAL] != MISSING)
                                            & (G[:, PATERNAL] != G[:, MATERNAL])))
                if unphased_hets[0].size:
                    if self.debug:
                        if count >= num_hom_haps:
                            print 'After hom, items on queue %d' % (q.qsize(),)
                            pass
                        print 'unphased_hets', unphased_hets
                        print 'si', si
                        print 'index i[unphased_hets]', np.where(i)[0][unphased_hets]
                        print 'H_unphased', H[i][unphased_hets]
                        print 'G_unphased', G[unphased_hets]
                    H_unphased = H[i][unphased_hets]
                    H_phased = H_unphased.copy()
                    complete_haplotype_partial(H_phased, G[unphased_hets])
#                    if self.debug:
#                        print 'Phasing hets'
#                        print 'i', np.where(i) 
#                        print 'H_unphased', H_unphased
#                        print 'G of unphased_hets', G[unphased_hets]
                    newly_phased_alleles = np.where(H_phased != H_unphased)[1]
                    self.h[si[unphased_hets]] = H_phased[:]
#                    if self.debug:
#                        print 'After phasing H_unphased', self.h[s[unphased_hets]]
#                        print 'unphased_hets', s[unphased_hets]
#                        print 'newly_phased_alleles', newly_phased_alleles
                    # Append the new data we can now make use of to the queue
                    if self.debug:
                        print 'After phasing them'
                        print 'H_phased  ', H_phased
                    for hap in zip(si[unphased_hets], newly_phased_alleles):
                        if self.debug:
                            print 'Adding phased het haplotypes to queue', hap
                        q.put(hap)
        self.__override_training_imputed_by_genotypes()
        if self.debug_sample >= 0:
            j = self.debug_sample
            if j in self.training_index: print 'Final g[%d] = %s, h[%d] = %s' % (self.training_index[j], repr(self.g[self.training_index[j]]), j, repr(self.h[j]))
            else: print 'h[%d] = %s' % (j, repr(self.h[j]))

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __override_training_imputed_by_genotypes(self):
        '''Override training set imputed genotypes by their original genotypes in order not to lose data.'''
        t = self.training_sample_index 
        G, H = self.g, self.h[t]
        need_override = np.where((H[:, 0] != MISSING) + (H[:, 1] != MISSING) < (G[:, 0] != MISSING) + (G[:, 1] != MISSING))[0]
        if need_override.size:
            # print 'Overriding training imputed by genotypes'
            if self.debug:
                for i in need_override:
                    print i, H[i], 'replaced with ', G[i]
            self.h[t[need_override]] = G[need_override]
            self.hap_type[t[need_override]] = UNPHASED

    def __phase_hom(self):
        '''Phase all homozygous training samples; place in corresponding locations in h = imputed_data.
        Set the training sample ID subsets hom, non_hom.'''
        hom_index = np.where(im.gt.is_homozygous(self.g)[:])
        hom = self.training_sample_index[hom_index]
        self.h[hom, :] = self.g[hom_index]
        return hom
