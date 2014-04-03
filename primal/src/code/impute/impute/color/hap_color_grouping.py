'''
============================================================
Old implementation of haplotype coloring by grouping
common hap segments = connected components of a graph.

Created on August 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import util, numpy as np, networkx as nx, itertools as it, impute as im, matplotlib.pyplot as P, operator

#---------------------------------------------
# Constants
#---------------------------------------------

#---------------------------------------------
# Methods
#---------------------------------------------
def to_group_to_color(segment_set, samples=None, segment_gap=25, snp_range=None):
    '''Group a SegmentSet object to segment sets of unique color.'''
    return group_to_color([(s.snp, s.samples) for s in segment_set], samples=samples, segment_gap=segment_gap, snp_range=snp_range)

def group_to_color(segments, samples=None, segment_gap=25, snp_range=None):
    '''Given a -list- of key-value pairs ((start, stop), (v1,v2)), where the key denotes
    a segment [start,stop) of the infinite integer lattice on which v1 and v2 are equal,
    
    Generate a haplotype coloring plot from a list of IBD segments. samples is the set of 
    haplotype identifiers. if pair_gap is specified, haplotypes pairs are separated by this
    amount [pixels], assuming that their number is even. segment_gap is the # of SNPs to drop
    from each segment''s margins for the purpose of defining colors. This allows a slight
    overlap between segments without making them the same color, which is usually undesirable.'''
    
    # Define colors using slightly-smaller portions of the original segments
    segment_groups = __group_to_color([(__dilute_segment(x[0], segment_gap), x[1]) for x in segments], samples, snp_range=snp_range)

    # Translate back result to the original segments. Use sub-segments, dangling nodes of the
    # original segments (note: there may be more irrelevant dangling nodes in the diluted segments)
    sub_segments, _, _, segment_to_local_values, dangling_local_values = \
    im.segment.form_intersections(segments, True, samples=samples, snp_range=snp_range)
    
    # Translate each connected component from a list of segments to a list of local values
    # (i.e., haplotype parts of the same color, for each color)
    return sub_segments, np.array([list(util.union_all(*(segment_to_local_values[x] for x in group))) 
                                   for group in segment_groups] + [[x] for x in dangling_local_values])

haps_of_color = lambda d, c: [(d[0][k], y[1]) for k, x in enumerate(d[1]) for y in x if y[0] == c]

#---------------------------------------------
# Plot Methods
#---------------------------------------------
def plot_hap_coloring(segment_set, haps=None, title=None, xlabel=None, ylabel=None, y=None,
                      pair_gap=0, linewidth=6, segment_gap=25, snp_range=None):
    '''Generate a haplotype coloring plot from a list of IBD segments. haps is the set of 
    haplotype identifiers. if pair_gap is specified, haplotypes pairs are separated by this
    amount [pixels], assuming that their number is even.
    
    segment_gap is the # of SNPs to ignore from each segment''s margins for the purpose of
    defining colors. This allows a slight overlap between segments without making them the
    same color, which is usually undesirable.'''
    # Generate haplotype colors 
    plot_hap_coloring_from_colors(segment_set.to_group_to_color(samples=haps, segment_gap=segment_gap, snp_range=snp_range),
                                  haps=haps, title=title, xlabel=xlabel, ylabel=ylabel, y=y, pair_gap=pair_gap,
                                  linewidth=linewidth, segment_gap=segment_gap, snp_range=snp_range)

def plot_hap_coloring_from_colors(d, haps=None, title=None, xlabel=None, ylabel=None, y=None,
                      pair_gap=0, linewidth=6, segment_gap=25, snp_range=None, colors=None):
    '''Generate a haplotype coloring plot from coloring scheme d.'''
    # Generate haplotype colors 
    sub_segments, groups = d
    haps = haps if haps is not None \
    else sorted(util.union_all(*(map(operator.itemgetter(1), x) for x in groups)))

    # Prepare axes
    colors = colors if colors is not None else im.plot.colors.get_colors()
    num_haps = len(haps)
    hap_index = dict((hap, k) for k, hap in enumerate(haps))
    custom_label = y is not None
    y = y if y is not None else map(repr, haps)
    xmin, xmax = sub_segments[0][0] - 1, sub_segments[-1][1] + 1
    x_scaled = lambda x: (1.0 * (x - xmin)) / (xmax - xmin)
    
    P.clf()
    P.hold(True)
    if title is not None:
        P.title(title)
    P.xlim([xmin, xmax])
    P.xlabel(xlabel if xlabel else 'SNP #')
    P.ylabel('Sample')
    hap_ticks = np.arange(0, num_haps)
    if pair_gap > 0:
        hap_ticks += np.array(list(it.chain.from_iterable(map(lambda x: (x, x), xrange(num_haps / 2)))))
    # hap_ticks = list(reversed(hap_ticks))
    P.yticks(hap_ticks, y)
    P.ylim([hap_ticks[0] - 0.5, hap_ticks[-1] + 0.5])
    
    for ((_, group), color) in it.izip(enumerate(groups), colors):
        # Draw lines for all members of the group using the same color
        for k, hap in group:
            print 'region %-2d [%-4d,%-4d]  x %.2f:%.2f  hap (%-4d, %d)%s  color (%.2f, %.2f, %.2f)' % \
            (k, sub_segments[k][0], sub_segments[k][1],
             x_scaled(sub_segments[k][0]), x_scaled(sub_segments[k][1]),
             hap[0], hap[1], ' ylabel %-10s' % (y[hap_index[hap]],) if custom_label else '',
             color[0], color[1], color[2])
            P.axhline(y=hap_ticks[hap_index[hap]],
                      xmin=x_scaled(sub_segments[k][0]),
                      xmax=x_scaled(sub_segments[k][1]),
                      linewidth=linewidth, color=color)
    P.show()
    return sub_segments, groups

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __dilute_segment(segment, segment_gap):
    '''Chop a margin of size segment_gap SNPs on either side of the segment tuple representing
    the SNP segment [a,b).'''
    start, stop = segment
    # Chop only if segment is not too small so that it would be empty after dilution
    return segment if stop - start <= 2 * segment_gap else (start + segment_gap, stop - segment_gap)

def __group_to_color(segments, samples=None, snp_range=None):
    '''Given a -list- of key-value pairs ((start, stop), (v1,v2)), where the key denotes
    a segment [start,stop) of the infinite integer lattice on which v1 and v2 are equal,
    Return a list of lists, each of contains segments of the same color (IBD sharing).'''
    (sub_segments, intersections, value_to_segments, _, _) = \
    im.segment.form_intersections(segments, True, samples=samples, snp_range=snp_range)
    # Build a graph G where the nodes = segments and edges = (segments intersect AND their sample
    # sets intersect). G's connected components are groups, where group is a set of segments of the
    # same color.
    return nx.connected_components(nx.from_edgelist(it.chain.from_iterable(it.product(sharing_segments, sharing_segments)
              for sharing_segments in (util.union_all(*(value_to_segments[i][x] for x in component))
                                       for i in xrange(len(sub_segments))
                                       for component in nx.connected_components(nx.Graph(intersections[i]))))))
    
####################################################################################
