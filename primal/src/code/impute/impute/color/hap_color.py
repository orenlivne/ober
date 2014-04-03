'''
============================================================
Old implementation of haplotype coloring by grouping
common hap segments = connected components of a graph.

Created on August 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, itertools as it, impute as im, matplotlib.pyplot as P
from impute.ibd.segment import START, STOP
from statutil import POSITIVE_TO_MINUS_1_TO_1

#---------------------------------------------
# Constants
#---------------------------------------------
'''All partitions of four paternal haplotypes into subsets of two.''' 
PATERNAL_HAP_PARTITIONS = map(lambda x: tuple(tuple(y) for y in x), frozenset([frozenset([frozenset(x), frozenset(np.arange(4)) - frozenset(x)]) for x in it.combinations(np.arange(4), 2)]))

'''Useful vectorized functions.'''
_VALUE_OR_ZERO = np.vectorize(lambda x, y: x if y else 0)
_IN_SET = np.vectorize(lambda x, y: x in y)

#---------------------------------------------
# Methods
#---------------------------------------------
def hap_colors(haps, segments, poo_phase=None, max_colors=100, max_sweeps=2):
    '''Color a list of haplotypes identified by ''haps'' that share the IBD segments in ''segments.''
    Use at most ''num_colors+1'' (num_colors+1 represents the ''unknown region'' color = doesn''t
    belong to any of the other num_colors colors) colors.'''
    # Initialize palette object, map haps and disjoint regions to ordinal integers  
    hap_label = dict((v, k) for k, v in enumerate(haps))    
    disjoint_segments = im.segment.SegmentSet(im.segment.Segment(x.snp, [hap_label[y] for y in x.samples], x.bp) for x in segments.to_group_to_disjoint())
    regions = sorted(set([x.snp for x in disjoint_segments]))
    region_lengths = np.array(map(lambda x: x[STOP] - x[START], sorted(set([x.bp for x in disjoint_segments]))), dtype=np.long)  # must match regions'' ordering and does in this case
    region_num = dict((v, k) for k, v in enumerate(regions))
    pa = Palette(haps, regions, region_lengths, max_colors + 1)
    (m, n), c = pa.shape, pa.color
    
    # Initialize each haplotype to its own color 
    for i in xrange(m): c[i] = i

    # Pass 1: collapse all hap#'s sharing an IBD segment to the minimum hap# among them
    # For each region, keep the set of IBD colors
    for x in disjoint_segments:
        region = region_num[x.snp]
        y = np.array(list(x.samples))
        c[y, region] = min(y)
        pa.cliques[region].append(frozenset(x.samples))
        print region, y, min(y), pa.cliques[region]
    
    # Identify primary colors. Recode colors in descending abundance (0=most abundant, ..., n-1=least abundant)
    colors = np.argsort(map(pa.color_coverage, xrange(pa.num_haps)))[::-1]
    K = len(colors)
    new_color = [0] * K
    for i in xrange(K): new_color[colors[i]] = i
    c[:] = np.vectorize(new_color.__getitem__)(c)
    # Recode clique color numbers
    for region in xrange(n):
        cliques = filter(lambda x: x, ([new_color[y] for y in x if new_color[y] < max_colors] for x in pa.cliques[region]))
        print 'Region', region
        if region == 50:
            pass
        print 'Only primary colors', cliques
        # Add lone-color cliques
        for k in xrange(max_colors):
            if not any(k in clique for clique in cliques): cliques.append([k])
        pa.cliques[region] = [np.array(clique) for clique in cliques]
        print 'With lone colors', pa.cliques[region]
        # Index cliques
        for index, clique in enumerate(cliques): pa.clique_index[clique, region] = index
    
    # Pass 2: try to minimize non-primary colors
    r, s = np.empty((K,), dtype=np.int), np.empty((max_colors,), dtype=np.int)
    for j in xrange(n - 1):
        # print 'Column', j
        p, q = c[:, j], c[:, j + 1]
        # r[k] = number of continuous rows j->j+1 of color k
        r.fill(0)
        for i in xrange(m):
            if p[i] == q[i]: r[p[i]] += 1
        for k in xrange(max_colors, K):
            # Determine best primary color b to swap non-primary color k with.
            # c minimized the number of discontinuities j->j+1 after the swap.
            where_k = (q == k)
            pk = p[where_k]
            s.fill(-1)
            for b in xrange(max_colors): s[b] = len(np.where(pk == b)[0]) - r[b]
            b = np.argmax(s)
            if s[b] > 0:
                # print 'Swapping', k, 'with', b
                # Found primary color, swap
                where_b = (q == b)
                c[where_k, j + 1] = b
                c[where_b, j + 1] = k
        # c = best primary color
    
    # Step 3: Replace all non-primary colors by a single color
    c[c >= max_colors] = max_colors
    # Replace grays by available primary colors whenever possible
    primary_colors = frozenset(xrange(pa.num_colors - 1))
    non_primary = pa.num_colors - 1
    for j in xrange(n):
        cj = c[:, j]
        # print j, cj
        available_colors = list(primary_colors - set(cj))
        # print 'available_colors', available_colors
        cj_non_primary = np.where(cj == non_primary)[0]
        # print 'cj_non_primary', cj_non_primary
        num_assignable = min(len(cj_non_primary), len(available_colors))
        if num_assignable > 0:
            # print 'Row', i, 'col', j, 'Assigning available', available_colors, 'to non_primary locations', cj_non_primary
            c[cj_non_primary[:num_assignable], j] = available_colors[:num_assignable]
            
    # Step 4: relaxation sweeps to minimize # recombinations    
    for _ in xrange(max_sweeps):
        pa_new = relax(pa)
        if np.all(pa_new == pa.color): break
        pa.color = pa_new  
    return pa

def relax(pa):
    '''Relaxation sweep to minimize #recombinations.'''
    c, J = pa.color.copy(), pa.num_regions
    for i in xrange(pa.num_haps):
        for j in xrange(J):
            possible_values = pa.clique_of(c[i, j], j) 
            if len(possible_values) > 1:
                # Set c[i,j] to the value x that minimize the sum of recombinations (discontinuities)
                # with fixed left and right neighbors on haplotype i (row i of the matrix c).
                diff_measure = np.zeros_like(possible_values)
                if j > 0: diff_measure += [int(c[i, j - 1] != x) for x in possible_values]
                if j < J - 1: diff_measure += [int(c[i, j + 1] != x) for x in possible_values]
                c[i, j] = min(zip(diff_measure, possible_values))[1]
    return c

#---------------------------------------------
# Plot Methods
#---------------------------------------------
def plot_palette(pa, title=None, xlabel=None, y=None, show_regions=False, show_primary_haps=False,
                 pair_gap=0, linewidth=6, snp_range=None, colors=None):
    '''Generate a haplotype coloring plot from coloring palette pa.'''
    # Generate haplotype colors 
    haps, regions, c = pa.haps, pa.regions, pa.color
    
    colors = colors if colors is not None else list(it.islice(im.plot.colors.get_colors(), pa.num_colors))
    # Last color: always gray
    colors[-1] = (0.4, 0.4, 0.4)

    # Prepare axes
    num_haps = len(haps)
    y = y if y is not None else map(repr, haps)
    xmin, xmax = regions[0][START], regions[-1][STOP]
    x_scaled = lambda x: (1.0 * (x - xmin)) / (xmax - xmin)
    
    P.clf()
    P.hold(True)
    if title is not None: P.title(title)
    P.xlim([xmin, xmax])
    P.xlabel(xlabel if xlabel else 'SNP #')
    P.ylabel('Sample')
    hap_ticks = np.arange(0, num_haps)
    if pair_gap > 0:
        hap_ticks += np.array(list(it.chain.from_iterable(map(lambda x: (x, x), xrange(num_haps / 2)))))
    P.yticks(hap_ticks, y)
    ymin, ymax = hap_ticks[0] - 0.5, hap_ticks[-1] + 0.5
    P.ylim([ymin, ymax])

    # Show region boundaries in dashed lines
    if show_regions:
        boundaries = sorted(set([y for x in pa.regions for y in x]))
        print boundaries
        for k, boundary in enumerate(boundaries[1:], 1):
            P.axvline(x=boundary, linestyle='dotted', linewidth=1, color='k')
            P.text(0.5 * (boundaries[k - 1] + boundaries[k]), hap_ticks[0] - 0.4, '%d' % (k - 1,))
            
    # Draw each segment in its corresponding color
    for i, j in it.product(xrange(pa.num_haps), xrange(pa.num_regions)):
        # Draw lines for all members of the group using the same color
        # print '[%d,%d] x %.2f:%.2f y=%.2f color %d' % (regions[j][START], regions[j][STOP], x_scaled(regions[j][START]), x_scaled(regions[j][STOP]), hap_ticks[i], c[i, j])
        P.axhspan(xmin=x_scaled(regions[j][START]), xmax=x_scaled(regions[j][STOP]),
                  ymin=hap_ticks[i] - 0.5 * linewidth, ymax=hap_ticks[i] + 0.5 * linewidth,
                  color=colors[c[i, j]])
    
    # Draw where primary haplotypes have been identified
    if show_primary_haps:
        primary = pa.color_sequence()
        for k, seq in enumerate(primary): 
            for i, j in seq:
                P.axhspan(xmin=x_scaled(regions[j][START]), xmax=x_scaled(regions[j][STOP]),
                          ymin=hap_ticks[i] + 0.9 * linewidth, ymax=hap_ticks[i] + 1.1 * linewidth,
                          color=colors[k])
    P.show()

####################################################################################
class Palette(object):
    '''Defines a color palette over haplotypes. Each color represents an ancestral haplotype
    with which segments on haplotypes are IBD.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, haps, regions, region_lengths, num_colors):
        '''Initialize a color palette for haplotypes labeled by the unique identifiers ''haps''.
        Colors are uniform within each segment in the segment set ''regions''. Use at most
        num_colors colors (num_colors-1 primary + one non-primary)'''
        self.haps = haps
        self.regions = regions
        self.region_lengths = region_lengths
        self.total_length = sum(region_lengths)
        # Matrix of color numbers in each haplotype and region
        self.num_haps = len(haps)
        self.num_regions = len(regions)
        self.num_colors = num_colors
        self.color = np.zeros((self.num_haps, self.num_regions), dtype=np.int)
        # Dictionary of IBD color groups (cliques) in each region
        self.cliques = [list() for _ in xrange(self.num_regions)]
        # Entry [i,j] is the index of the clique in self.cliques[j] in which i lies
        self.clique_index = np.zeros((self.num_colors, self.num_regions), dtype=np.int)
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        s = 'Palette[#haps=%d, #regions=%d, #colors=%d]\n' % (self.num_haps, self.num_regions, self.num_colors)
        s += 'Haps %s\n' % (repr(self.haps),)
        s += 'SNP Regions %s\n' % (repr(self.regions),)
        s += 'Region lengths %s\n' % (repr(self.region_lengths),)
        s += 'Cliques:\n' + repr(self.cliques)
        s += '\nClique index:\n' + repr(self.clique_index)
        s += '\nColors:\n' + repr(self.color)
        return s
    
    @property
    def shape(self): return self.num_haps, self.num_regions
    
    def color_coverage(self, color):
        '''Return the total length of genomic regions covered by color number ''color''.'''
        return sum(self.region_lengths[np.where(self.color == color)[1]])

    def clique_of(self, k, j):
        '''Return the clique (list of IBD colors) of color k, region j.'''
        return self.cliques[j][self.clique_index[k, j]]

    def color_sequence(self):
        '''Return sequence of hap and region coordinates of each ''primary-colored'' haplotype.'''
        c = self.color
        K = self.num_colors - 1
        seq = [list() for _ in xrange(K)]
        for j in xrange(self.num_regions):
            # print 'col', j
            found_color = np.zeros((K,), np.bool)
            for i in xrange(self.num_haps):
                entry = c[i, j]
                if entry < K and not found_color[entry]:
                    found_color[entry] = True
                    seq[entry].append((i, j))
            # If hap not found but is IBD with a found color, append to its sequence
            # print 'found_color', found_color
            for k in np.where(~found_color)[0]:
                try:
                    l = (y for y in self.clique_of(k, j) if found_color[y]).next()
                    # print 'Can assign', k, 'to', l
                    found_color[k] = True
                    seq[k].append(seq[l][-1])
                except StopIteration: pass  # No found IBD color exists
        return seq

    def color_sequence_coverage(self, color_nums):
        '''Calculate the % of each haplotype in haps covered by the color sequence numbers
        (=paternal haplotype numbers) in the array ''color_nums''.'''
        # t = array of paternal haplotype color sequences
        # ibd_colors[i] = all colors that are IBD with ''colors'' in region[i] 
        t, ibd_colors = self.color_sequence(), [set() for _ in xrange(self.num_regions)]
        for c, j in ((c, j) for c in color_nums for _, j in t[c]): ibd_colors[j] |= set(self.clique_of(c, j))
        lengths, total_length = np.tile(self.region_lengths, (self.num_haps, 1)), sum(self.region_lengths)
        return _VALUE_OR_ZERO(lengths, _IN_SET(self.color, np.tile(ibd_colors, (self.num_haps, 1)))).sum(axis=1).astype(float) / total_length
        
    def parent_coverage_bp(self, color):
        '''Return the total base pairs of parent haplotype number ''color'' that can be recovered.'''
        return sum(self.region_lengths[region] for region in np.unique(np.where(self.color == color)[1]))
        #return [sum(self.region_lengths[v] for _, v in x) / float(self.total_length) for x in self.color_sequence()]
        
def coverage_ratio(pa, color_nums1, color_nums2):
    '''Given two groups of paternal haplotype color numbers color_nums1, color_nums2, calculate
    the ratio of coverage of each child chromosome pair, assuming the child POO phase is 1 and -1.
    Note: child haplotypes are assumed to be contiguous in the pa.haps array, i.e., child k = 
    haps[2*k],hap[2*k+1].'''
    a, b = pa.color_sequence_coverage(color_nums1), pa.color_sequence_coverage(color_nums2)
    a0, a1, b0, b1 = a[::2], a[1::2], b[::2], b[1::2] 
    return np.array([(a0 + b1) / (a1 + b0 + 1e-16), (a1 + b0) / (a0 + b1 + 1e-16)])

def best_hap_alignment_to_colors(pa):
    '''Align child POO phases by considering all possible partitions of paternal haplotypes.
    Returns an array indicating whether to flip a child or not, and the separation measure between
    the two phases (>~0.2: phase=1, no need to flip haps; <~-0.2: phase=-1, flip haps).'''
    max_separation = 1e-16
    for paternal_colors, maternal_colors in PATERNAL_HAP_PARTITIONS:
        a = coverage_ratio(pa, paternal_colors, maternal_colors)
        separation = min(a.max(axis=0) / (a.min(axis=0) + 1e-16))
        if separation > max_separation: max_a, max_separation = a, separation
    return POSITIVE_TO_MINUS_1_TO_1(max_a[0] / (1e-16 + max_a[1])), POSITIVE_TO_MINUS_1_TO_1(max_separation), paternal_colors, maternal_colors
