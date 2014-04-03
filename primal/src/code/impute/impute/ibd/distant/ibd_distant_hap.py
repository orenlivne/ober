'''
============================================================
IBD segment identification between two haplotypes.

Created on January 25, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, impute as im, itertools as it, itemutil, time
from impute.ibd.segment import START, STOP, Segment
from impute.ibd.distant.ibd_hmm_hap import prob_ibd_hmm_from_raw_input, IbdProblem
from multiprocessing import Pool, Manager
#from util import write_with_lock

#---------------------------------------------
# Methods
#---------------------------------------------
def sample_segments((problem, i, j, params, lock)):
    '''Return IBD segments between phased samples i,j.'''
    # write_with_lock('IBD Segments between %d, %d' % (i, j), lock)
    segments = sum((hap_segments(IbdProblem(problem, (i, a), (j, b), None, params))
                    for a, b in it.product(im.constants.ALLELES, im.constants.ALLELES)),
                    im.segment.SegmentSet([]))
    # write_with_lock('\tfound %d segments\n' % (segments.length,), lock)
    return segments

def among_samples_segments(problem, s, params):
    '''Return IBD segments among samples in the set s. Only supports a serial run for now.'''
    return sum((sample_segments((problem, i, j, params, None)) 
               for i, j in it.product(s, s) if i < j), im.segment.SegmentSet([]))

def between_samples_segments(problem, s, t, params, num_processes=1):
    '''Return IBD segments between samples in the set s and samples in the set t.
    s,t should be disjoint (if not, duplicate pairs will be present in the result set).
    If num_processes > 1, uses a multiprocessing pool with num_processes parallel processes.'''
    pairs = list(it.product(s, t))
    if num_processes == 1:  # Serial
        return sum((sample_segments((problem, i, j, params, None)) for i, j in pairs), im.segment.SegmentSet([]))
    else:  # Parallel
        manager = Manager()
        lock = manager.Lock()
        
        # Map phase
        po = Pool(processes=num_processes)
        start = time.time()
        res = po.map(__hap_segments_task, ((IbdProblem(problem, (i, a), (j, b), None, params), lock)
                     for i, j in pairs for a, b in it.product(im.constants.ALLELES, im.constants.ALLELES)))
        t = time.time() - start
        
        # Reduce phase
        total_segments = sum(res, im.segment.SegmentSet([]))
        if params.debug:
            print 'Total %d hap pairs, %d segments' % (4 * len(pairs), total_segments.length)
            print 'Elapsed Time: %.3f sec (%.3f sec/pair)' % (t, t / (4 * len(pairs)))
        return total_segments

'''Return IBD segments between pairs of samples in the tuple list s.'''
paired_samples_segments = lambda problem, s, params: \
sum((sample_segments((problem, i, j, params, None)) for (i, j) in s), im.segment.SegmentSet([]))

'''A utility method that conveniently accepts a problem object directly. Useful in a non-parallel
environment where we don''t need to serialize method arguments.'''
hap_segments_from_problem = lambda problem, hap1, hap2, params: \
hap_segments(IbdProblem(problem, hap1, hap2, None, params))

def hap_segments(h):
    '''Return IBD segments between two phased sample haplotypes (i,a),(j,b) encapsulated by
    the IbdProblem object ''h''.'''
    # Load data
    (i, a), (j, b), debug, num_snps, frame = h.hap1, h.hap2, h.params.debug, h.num_snps, h.snps
    if debug:
        print '-' * 70
        print 'IBD Segments between (%d,%d), (%d,%d)' % (i, a, j, b)
        print '-' * 70
        
    # Calculate IBD posterior using the largest frame
    prob_ibd = prob_ibd_hmm_from_raw_input(h)
    is_ibd = prob_ibd > 0.0
    if debug:
        # print 'IBD mask', is_ibd.astype(np.uint)
        print 'Frame IBD segments', [(frame[x[START]], frame[x[STOP]] if x[STOP] < len(frame) else num_snps)
                                     for x in im.segment.segments_with_value(is_ibd, True, 3)] 
    
    # Output sufficiently-long segments
    min_len_sub_segment = h.params.min_len_sub_segment  # If a segment is identified by HMM and is above params.min_len, allow smaller sub-segments in which the haps fit @UnusedVariable
    bp, cm, len_unit = h.bp, h.cm, h.params.len_unit
    coord = cm if (len_unit == 'cm') else bp / (1.0 * im.constants.MEGA_BASE_PAIR) 
    pair = [h.hap1, h.hap2]
    long_segments = im.segment.SegmentSet([])
    coord_bp = lambda x, y: (bp[x], im.segment.stop_bp(bp, y, num_snps))
    coord_cm = lambda x, y: (cm[x], im.segment.stop_bp(cm, y, num_snps))
    coord_unit = lambda x, y: (coord[x], im.segment.stop_bp(coord, y, num_snps))
    for fr, full_segment in ((fr, __frame_segment_to_snp(frame, fr, num_snps))
                             for fr in im.segment.segments_with_value(is_ibd, True)):
        # Interpolate prob_ibd from frame SNPs to all SNPs in the segments
        start, stop = full_segment[START], full_segment[STOP]  # Original indices of segment boundaries
        cm_start, cm_stop = coord_unit(start, stop)
        if debug:
            print '%d:%d coordinates %f:%f' % (start, stop, cm_start, cm_stop)
        if cm_stop - cm_start >= h.params.min_len:  # Threshold should be really small to catch everything, since we believe the IBD HMM to give the correct answer
            fr_start, fr_stop = fr[START], fr[STOP]  # Frame indices of segment boundaries
            confidence = np.interp(coord[start:stop], coord[frame[fr_start:fr_stop]], prob_ibd[fr_start:fr_stop])
            # Add entire segment
            # long_segments += Segment((start, stop), pair, (bp_start, bp_stop), confidence=confidence)
            # Check if haplotypes fit in ALL SNPs, not just SNP frames. Restrict to largest
            # sub-segment of consecutive "1"s found in [start,stop] (allowing genotype error mismatches)
            d = h.d[start:stop]
            d_filtered = im.ibd.filter_diff(d, 'median', 10)
            if debug:
                np.set_printoptions(threshold=np.nan)
                print '    %-4s %-4s %-4s %-4s %-4s' % ('SNP', 'h1', 'h2', 'd', 'filt')
                print np.concatenate((np.arange(start, stop)[np.newaxis].transpose(),
                                      h.all_h1[start:stop][np.newaxis].transpose(),
                                      h.all_h2[start:stop][np.newaxis].transpose(),
                                      d.astype(np.int)[np.newaxis].transpose(),
                                      d_filtered.astype(np.int)[np.newaxis].transpose()
                                      ), axis=1)

            # Use all sufficiently-long sub-segments
            ind = [v for k, v in itemutil.groupby_with_range(d_filtered) if k]
            s = map(lambda x: x + start, ind)  # segments of ones
            sub_segments = [segment for segment in
                            (Segment((int(s[i][START]), int(s[i][STOP])), pair, coord_bp(s[i][START], s[i][STOP]),
                                     confidence=confidence[ind[i][START]:ind[i][STOP]],
                                     cm=coord_cm(s[i][START], s[i][STOP]))
                             for i in xrange(len(ind)) if s[i][STOP] > s[i][START])
                             if (segment.length_cm if (len_unit == 'cm') else segment.length) >= min_len_sub_segment]
            if stop > start: long_segments += sub_segments
            if debug:
                print 'segments of ones', s
                print 'ind', ind
                print 'sub-segments', sub_segments
        else:
            if debug: print 'Segment shorter than %.2f threhold' % (h.params.min_len,)
    if debug:
        print 'Long-enough segments (segment >= %.2f %s, sub-segment >= %.2f %s):' % \
        (h.params.min_len, len_unit, h.params.min_len_sub_segment, len_unit)
        print long_segments.pprint_segments(show_bp=True) if long_segments else '\t-'
        print ''
    return long_segments

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __frame_segment_to_snp(frame, segment, num_snps):
    '''Convert a segment=[a,b) in frame coordinates to SNP indices of the original problem.'''
    return (frame[segment[START]], frame[segment[STOP]] if segment[STOP] < len(frame) else num_snps)

def __hap_segments_task((h, lock)):
    '''Same as hap_segment - a convenient wrapper for a parallel run. Accepts an IbdProblem and
    a processor pool lock.'''
    # write_with_lock('IBD Segments between (%d,%d), (%d,%d)\n' % (h.hap1[0], h.hap1[1], h.hap2[0], h.hap2[1]), lock)
    return hap_segments(h)
