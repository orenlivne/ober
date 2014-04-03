#!/usr/bin/env python
'''
============================================================
Impute SNPs from 98 CGI WGS samples to all Hutt samples.

Created on November 30, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, os
from impute.data.constants import MISSING, NUM_CHROMOSOMES
from impute.tools.genotype_tools import empty_errors_array
 
def merge_errors(errors1, errors2):
    '''Merge two error arrays. TODO: replace by a += operator.'''
    return np.concatenate((errors1, errors2), axis=1)

def phase_at_snp_using_training_set(snp_index, ibd, h, t, training_sample_set):
    '''Phase using IBD segments containing the training_sample_set.'''
    errors = empty_errors_array()
    bp = t.snp['base_pair'][snp_index]
    segments = im.imputation.ibd_lookup.segments_ibd_at_bp(ibd, (bp, bp), samples=training_sample_set)
    for segment in segments:
        s = np.array(list(segment.samples))
        hap_ids = np.array([x[0] for x in segment.samples])
        i = np.in1d(hap_ids, training_sample).nonzero()[0]
        (i1, i2) = s[i].transpose()
        alleles = h[snp_index, i1, i2]
        known = alleles.nonzero()[0]
        if known.size:
            # There exist phased training samples
            hap_ids = hap_ids[i[known]]
            alleles = alleles[known]
            if np.diff(alleles).nonzero()[0]:
                # Contradicting alleles, flag errors
                print 'Error: contradicting alleles at snp index', snp_index, 'hap_ids', hap_ids, 'alleles', alleles
                e = np.zeros((2, len(hap_ids)), dtype=np.uint)
                e[0, :] = snp_index
                e[1, :] = hap_ids                
                errors = np.concatenate((errors, e), axis=1)
            elif alleles.size:
                # Copy known alleles to all IBD haplotypes' h-entries
                (i1, i2) = s.transpose()
                print 'phasing at i1', i1, 'i2', i2, 'allele value', alleles[0]
                h[snp_index, i1, i2] = alleles[0]
    return errors

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Load IBD dictionary from final phasing result (merge all chromosomes and parts to one big dictionary).
    # This is not currently only a pairwise IBD dictionary, ot a a global IBD 
    phasing_dir = os.environ['OBER'] + '/out/hutt/phasing.20121130'
    p = np.load('%s/reduce/result/hutt_phased.stats.npz' % (phasing_dir,))
    pedigree = p['pedigree'][0]
    info = p['info'][0]
    ibd = [im.ibd.segment.SegmentSet(sum((x.ibd for x in info[chrom - 1]), [])) for chrom in xrange(NUM_CHROMOSOMES)]
    
    # Load training set genotypes
    t = im.io_genotype.read('npz', 'genotype', file=os.environ['OBER'] + '/data/impute/rare/rare.npz')
    tg = t.data 
    training_sample = np.array([pedigree.node_of[x] for x in t.sample_id])
    training_sample_set = set(training_sample)
    
    # Allocate imputed data structures
    imputed = im.factory.GenotypeFactory.new_instance('haplotype',
                                                      np.zeros((t.num_snps, pedigree.num_genotyped, 2), dtype=np.byte),
                                                      t.snp, t.sample_id)
    h = imputed.data
    
    # Phase all homozygous trainees; place in corresponding locations in h
    hom = np.where(im.gt.is_homozygous(tg)[:, :])
    h[hom[0], training_sample[hom[1]], :] = tg[hom]
    
    # Phase training set. For each SNP:
    # - Find all segments intersecting the SNP
    # - Construct global IBD segments-sets. This is currently a wasteful implementation. TODO: replace ibd by
    #   the global IBD dictionary to speed this part up
    # - Copy known alleles to 
    # - If known alleles differ, mark errors (at all alleles for now; TODO: use instead majority vote to find a single error among multiple alleles)
    errors = empty_errors_array()
    for snp_index in xrange(t.num_snps):
        chrom = t.snp['chrom'][snp_index]
        print '====== SNP %d: chr%d:%d, %s ======' % (snp_index, chrom, t.snp['base_pair'][snp_index],
                                                      t.snp['name'][snp_index])
        ibd_chrom = ibd[chrom - 1]
        # Phase using homozygous training samples 
        errors = merge_errors(errors, phase_at_snp_using_training_set(snp_index, ibd_chrom, h, t, training_sample_set))
        # Phase using het training samples that now have one allele determined          
        hh = h[snp_index, training_sample, :]
        het_to_phase = training_sample[np.where(((hh[:, 0] == MISSING) ^ (hh[:, 1] == MISSING)) & 
                                                (tg[snp_index, :, 0] != MISSING) & (tg[snp_index, :, 1] != MISSING))[0]]
        errors = merge_errors(errors, phase_at_snp_using_training_set(snp_index, ibd_chrom, h, t, het_to_phase))
