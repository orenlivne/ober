'''
============================================================
General-purpose algorithms related to Identity-By-Descent
(IBD) segment identification among haplotypes. Builds upon
the segment module.

Created on July 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, util
from scipy import ndimage
from impute.data import constants
from impute.data.constants import INDETERMINATE
from impute.tools import genotype_tools as gt, recode
from impute.ibd import diff
from impute.ibd.segment import Segment, SegmentSet
from impute.ibd import segment
from impute.tools.param import PhaseParam

#---------------------------------------------
# Methods
#---------------------------------------------
def filter_diff(d, error_filter, error_filter_length):
    '''Filter the hap difference array d using the filter of type error_filter.'''
    # d is a binary image. Filter genotype errors, manifested as spikes in d
    if error_filter == 'median':
        return ndimage.median_filter(d, error_filter_length)
        # filtered_diff = ndimage.median_filter(d_snps, error_filter_length, mode='constant', cval=0)
    elif error_filter == 'fast_max':
        center_value = 2 * error_filter_length + 1
        window = [1] * center_value
        window[error_filter_length] = center_value
        f = ndimage.convolve(d, np.array(window, dtype=int))
        return ((f > 1) & (f < center_value)).astype('int')
    elif error_filter == 'max_except_center':
        return util.max_except_center_filter(d, error_filter_length)
    else:
        raise ValueError('Unsupported error filter ''%s''' % (error_filter,))

def genotype_ibs_segments(genotype, id1, id2, snps,
                          error_filter='median', error_filter_length=5, margin=0.0,
                          min_ibs_len_snp=400, debug=False):
    '''Return Identical-by-State (IBS >= 1) segments between two genoypes of samples id1 and id2
    in the SNP range [snp[0],snp[1]) (if snp is a tuple) or the subset of SNPs, if snps is an array.
    
    See ibs_segments() for a description of optional parameters.'''
    num_snps = genotype.num_snps
    g = genotype.data
    g1 = recode.recode_single_genotype(g[snps, id1, :])
    g2 = recode.recode_single_genotype(g[snps, id2, :])
    d = (recode.ibs_state(g1, g2) == 0).astype(np.byte)

    # Consider informative or the specified SNPs only
    filtered_diff = filter_diff(d, error_filter, error_filter_length)
    error_snps = snps[np.nonzero(d - filtered_diff)[0]]
    
    # Detect edges as non-zero gradient points; output sufficiently long segments
    if np.size(filtered_diff) == 0:
        # No data to consider ==> no IBD intervals can be identified
        segments = []
    else:
        # Convert recombination locations to segments of no recombination; filter short segments
        bp = genotype.snp['base_pair']
        #print segment.segments_with_value(filtered_diff, 0, min_ibs_len_snp)
        segments = [Segment(((x[0], x[1])), [id1, id2], (bp[x[0]], segment.stop_bp(bp, x[1], num_snps)),
                            error_snps=segment.in_segment(error_snps, x), collapse_to_set=False)
                    for x in segment.segments_with_value(filtered_diff, 0, min_ibs_len_snp)]
    
    # Cut segment margins
    if margin >= constants.SMALL_FLOAT:
        segments = [s for s in (s.middle_part(genotype.nearest_snp, bp, margin, collapse_to_set=False)
                                for s in segments) if s]

    # Restrict errors to those inside segments
    segment_set = SegmentSet(segments,
                             np.array(util.flattened_meshgrid(reduce(list.__add__, (s.error_snps.tolist() for s in segments)),
                                                              np.array([id1, id2])), dtype=int) \
                             if segments else gt.empty_errors_array())
    if debug:
        print 'ibs_segments()', segment_set
        print 'errors', segment_set.errors
    return segment_set

def ibs_segments(haplotype, id1, id2, hap1_type, hap2_type, snps=None, include_alt_phase=False,
                 error_filter='median', error_filter_length=5,
                 length_bound=None, min_segment_length=INDETERMINATE, margin=0.0, debug=False):
    '''Return 1) Identical-by-State (IBS) segments separated by recombination events between two
    sample haplotypes (id1, hap1_type) and (id2, hap2_type). The 2-D output array's ith row format is
    
    (segment_start, segment_stop),
    (id1, hap1), (id2, hap2), 
    (segment_start_bp, segment_stop_bp, segment_length_in_bp, num_errors_in_segment) 
    
    The SNP range is [segment_start, segment_stop) where start=inclusive and stop is exclusive.
    2) List of het_snp indices at which there are likely genotype errors.
        
    Options:
    snps - list of SNPs to base the comparison on. For parent-child comparisons, these should
    be heterozygous SNPs in the parent's genotype, distinguishing its haplotypes
    and used to locate segments. For unphased-phased individuals, these should be the list of
    homozygous SNPs at the unphased individual (those that have data).
    If not specified, all SNPs are used.
    
    length_bound - minimum segment length bound type:
        None: no lower bound enforced 
        'base_pair': output segments of at least min_segment_length [base pair]
        'snp': output segments of at least min_segment_length consecutive SNPs out of the snps list.
               This is useful only if snps includes all SNPs (or is None) 
        *NOTE*: min_segment_length''s units are interpreted differently depending on length_bound.
         
    margin = fraction of segment to discard near the endpoints (margin/2 is removed from each side).'''
    
    if debug:
        print 'Computing IBD segments between haplotypes (%d,%d), (%d,%d); filter %s length %d' % \
        (id1, hap1_type, id2, hap2_type, error_filter, error_filter_length)
    d = diff.all_diffs(haplotype.data, id1, id2, hap1_type=hap1_type, hap2_type=hap2_type)[0]
    # Segment length, as defined by the input parameters 
    segment_length = lambda f: np.inf if not length_bound else (f.length if length_bound == 'base_pair' else f.num_snps)  # @UnusedVariable
    
    # Consider informative or the specified SNPs only
    snps = snps if snps is not None else haplotype.snp_range
    snps = np.intersect1d(snps, np.where(d != INDETERMINATE)[0])
    d_snps = d[snps]
    filtered_diff = filter_diff(d_snps, error_filter, error_filter_length)    
    error_snps = snps[np.nonzero(d_snps - filtered_diff)[0]]
    
    # Detect edges as non-zero gradient points; output sufficiently long segments
    bp = haplotype.snp['base_pair']
    num_snps = haplotype.num_snps
    if np.size(filtered_diff) == 0:
        # No data to consider ==> no IBD intervals can be identified
        segments = []
    else:
        deriv = ndimage.convolve(filtered_diff, [1, -1])    
        edge = np.where(deriv != 0)[0]
        initial_phase = hap1_type if filtered_diff[0] == 0 else 1 - hap1_type
        if debug:
            print 'initial_phase', initial_phase  # , 'edge', edge
        # Convert recombination locations to segments of no recombination; filter short segments
        segments = [f for f in (Segment(((x[0], x[1])), set(((id1, x[2]), (id2, hap2_type))),
                                        (bp[x[0]], segment.stop_bp(bp, x[1], num_snps)),
                                        error_snps=segment.in_segment(error_snps, x))
                                for x in segment.edges_to_segments(snps, edge, initial_phase,
                                                                   haplotype.num_snps, hap1_type))
                    if segment_length(f) >= min_segment_length]
    
    # Cut segment margins
    if margin >= constants.SMALL_FLOAT:
        segments = [s for s in (s.middle_part(haplotype.nearest_snp, bp, margin) for s in segments) if s]
    
    # Restrict errors to those inside segments
    segment_set = SegmentSet(segments,
                             np.array(util.flattened_meshgrid(reduce(list.__add__, (s.error_snps.tolist() for s in segments)),
                                                                        np.array([id1, id2])), dtype=int) \
                             if segments else gt.empty_errors_array())
    if debug:
        print 'ibs_segments()', segment_set
        print 'errors', segment_set.errors
    return segment_set

def concatenate_segments(segments):
    '''Concatenate segments and errors from a generator expression.'''
    (all_segments, all_errors) = ([], gt.empty_errors_array())
    for s in segments:
        all_segments += s
        all_errors = gt.concatenate_errors(all_errors, s.errors)
    return SegmentSet(all_segments, all_errors)

def phase_by_ibd(request, segments, consensus, phase=True, num_sweeps=2):
    # remove_autozygous=False, 
    '''A a template method to phase by computing IBD segments.''' 
    (debug, problem) = (request.params.debug, request.problem)
    
    segments = concatenate_segments(segments)
    
    if debug:
        print 'phase_by_ibd()'
        print 'Raw segments'
        print segments
        print 'IBD dictionary length before adding segments', problem.info.ibd.length
        print '#segments to add', segments.length
    # Save segments in IBD dictionary for future reference, e.g., imputation
    problem.info.ibd += segments
    if debug:
        print 'IBD dictionary length after adding segments', problem.info.ibd.length

    segments.group_to_disjoint()
#    if remove_autozygous:
#        segments.remove_duplicate_samples()
    if debug:
        print 'Grouped to disjoint sub-segments:'
        print segments
        print ''
    
    # Flag genotype errors found
    errors = segments.errors
    problem.genotype_error(errors[0], errors[1], 'in IBD segment calculation')
    if debug:
        print 'Raw errors %d' % (len(errors[0]),)
        
    if phase:
        # Impute in each segment: if a non-missing value is found, it is copied to all other
        # IBD-sharing individuals. This is cleanly implemented using max(h over samples).
        # Each segment may reveal additional genotype errors, if an individual is inconsistent
        # with an overwhelming majority vote.
        #
        # We sweep through all intervals twice to ensure full information propagation. As in
        # relaxation, the phasing information is propagated to neighbors on the IBD graph.
        # Thus, if say {(0,0),(1,0),(2,0)} and {(0,1),(3,0),(4,0)} are IBD hap sets over a
        # certain SNP segment, where 3 is homozygous and the rest are not and if the 0-1-2 set
        # is happened to be processed first, then 3 will propagate its allele to 0 and 4 only.
        # 0, however, has a full haplotype now, inferred from the (0,1) hap value and 0's genotype.
        # This value can now be used to phase 1 and 2 in the set 0-1-2. This will be caught by
        # the second iteration over segments.
        #
        # Note: the two-sweep loop can be processed in parallel for disjoint SNP segments, even
        # though it is not implemented that way at present.
        for sweep in xrange(num_sweeps):
            if debug:
                print '-' * 70
                print 'Sweep', sweep
                print '-' * 70
            for segment in segments:
                snp = segment.snp
                haps = segment.samples
                ordered_haps = np.array(list(haps))
                if debug:
                    print 'Phasing in IBD segment', snp
                    print ordered_haps
                _phase_in_ibd_segment(request.problem, snp, ordered_haps, consensus, debug=debug,
                                      params=request.params)

#---------------------------------------------
# Private Methods
#---------------------------------------------\
def _phase_in_ibd_segment(problem, snp, haps, consensus, debug=False, params=PhaseParam()):
    '''Impute haplotypes in an IBD-sharing sample set haps over a SNP array snp, or a SNP array
    segment [snp[0],snp[1]], if it is a tuple.
    
    Use the concensus functor to calculator a consensus haplotype and copy it to all other
    IBD-sharing samples. haps is a list of (sample,hap) tuples, where sample=sample ID and
    hap=allele (paternal/maternal) that are assumed to be IBD in this segment. 
    The other allele in each haplotype is inferred from the genotype.
    
    concensus = 'max': if a non-zero value is found, it is the concensus. This should be applied
    when all haps are certain, so all non-missing values should agree.
    This is cleanly implemented using max(h over samples). 
    
    concensus = 'majority': majority vote of non-missing values.'''
    
    # print 'Phasing in segment (%d,%d)' % (start, stop, )
    # problem, params = request.problem, request.params
    # snp_test_index = params.snp_test_index
    snps = np.arange(snp[0], snp[1]) if isinstance(snp, tuple) else snp
    common, hh = _compute_conensus(problem.h, snps, haps, consensus)

    # Flag samples that are inconsistent with the concensus as errors if we have enough haps
    # to support the evidence for errors. If there are too many errors, this is a dubious segment
    errors = _find_errors(problem, snps, haps, consensus, common, hh, params.min_consensus_samples)
    if debug: print 'Consensus errors (%d)' % len(errors[0]), errors
        # print 'Consensus errors %d' % len(errors[0])
    problem.genotype_error(errors[0], errors[1], 'IBD majority vote inconsistency')
    
    # Phase all haps: copy consensus haplotype to missing haplotype entries
    _phase_by_consensus(problem, snps, haps, common)

def _compute_conensus(h, snps, haps, consensus):
    hh = np.array([h[snps, sample, hap] for (sample, hap) in haps])
    consensus_func = _consensus_func(consensus)
    common = consensus_func(hh)
    return common, hh

def _consensus_func(consensus):
    if consensus == 'max':
        return _concensus_max
    elif consensus == 'majority':
        return _concensus_majority
    else:
        raise ValueError('Unsupported concensus type ''%s''' % (consensus,))

def _phase_by_consensus(problem, snps, haps, common):
    # Impute according to the common haplotype. Must be done after errors are filtered
    # so that we do not complete the haplotype, then flag an error in the first allele, but
    # leave behind a set second allele phased according to that first allele.
    g, h = problem.data  
    for sample, hap in haps:
        h[snps, sample, hap] = common
        gt.complete_haplotype_sample(h, g[snps, sample, :], snps, sample,
                                     error_status=problem.error[snps, sample])
#        gt.complete_haplotype(h[start:stop, sample, :], g[snps, sample, :],
#                              error_status=problem.error[snps, sample])

#    if debug and snp_test_index and start <= snp_test_index and snp_test_index <= stop:
#        print 'consensus', common[snp_test_index - start]
#        print 'Haps at SNP %d after phasing:' % (snp_test_index,)
#        ind = np.array(haps)
#        print h[max(start, snp_test_index - 5):min(stop, snp_test_index + 5), ind[:, 0], ind[:, 1]] 

def _get_current_error_status(problem, haps, snps):
    #return problem.error[snps][:, haps[:, 0]].transpose()
    # Much faster using http://stackoverflow.com/questions/14386822/fast-numpy-slicing
    return problem.error[snps[:, None], haps[:, 0]].transpose()

def _find_new_errors(common, hh, e):
        # Flag errors only if does not agree with consensus, is not missing and is not already an error
    return np.where((np.tile(common, (hh.shape[0], 1)) != hh) & (hh != constants.MISSING) & (e == constants.OK))

def _find_errors(problem, snps, haps, consensus, common, hh, min_consensus_samples):
    e = _get_current_error_status(problem, haps, snps)
#    if debug and snp_test_index and start <= snp_test_index and snp_test_index <= stop:
#        print 'Haps at SNP %d before phasing:' % (snp_test_index,)
#        np.set_printoptions(threshold=np.nan)
#        ind = np.array(haps)
#        print h[max(start, snp_test_index - 5):min(stop, snp_test_index + 5), ind[:, 0], ind[:, 1]] 
    if consensus == 'majority' and hh.shape[0] >= min_consensus_samples:
        # errors = np.where(np.tile(common, (num_haps,1)) != hh)
        # Convert to original coordinates: (SNP, sample)
        errors = _find_new_errors(common, hh, e)
        errors = (snps[errors[1]], haps[errors[0], 0]) if errors[0].size else errors
    else:
        errors = gt.empty_errors_array()
    return errors

def _concensus_max(h):
    '''If a non-zero value is found, it is the concensus. This should be applied
    when all haps are certain, so all non-missing values should agree.
    This is cleanly implemented using max(h over samples).''' 
    return np.max(h, axis=0)

def _concensus_majority(h):
    '''Majority voting of non-missing values. Not sure if this is the fastest vectorized implementation.''' 
    # - Make sure not to divide by 0 if all values are missing
    # - If the same number of votes is obtained for 1 and 2 alleles, round(2*(avg non-zero h)) will be 
    #   an odd number (e.g., equal votes of 1's and 2's ==> this is 3). Zero out these values so that
    #   the h-values will be different than the consensus and subsequently flagged as errors.  
    double_value = np.round(2 * np.sum(h, axis=0, dtype=float) / np.maximum(1e-15, np.sum(h != 0, axis=0))).astype(int)
    double_value[np.mod(double_value, 2) == 1] = constants.MISSING
    return double_value / 2

def _round_position(haplotype, bp_position):
    '''Round a base-pair position to the nearest snp and corresponding bp.'''
    snp = haplotype.nearest_snp(bp_position)
    bp = haplotype.snp['base_pair'][bp_position]
    return (snp, bp)
