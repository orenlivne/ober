'''
============================================================
Miscellaneous plot functions.

Created on August 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, os, util, impute as im, db_gene, matplotlib.pylab as P
from impute.ibd import diff
from impute.tools import genotype_tools as gt
from impute.data.constants import INDETERMINATE, PATERNAL, ALLELES, SMALL_FLOAT, NUM_CHROMOSOMES, ALLELE_LABELS
from impute.data import constants
from impute.data.constants import SNP, SAMPLE
from numpy.lib.type_check import isreal
from scipy import ndimage
from numpy.core.function_base import linspace

#---------------------------------------------
# Constants
#---------------------------------------------

# Haplotype color scheme 
# DEFAULT_COLORS = {0: 'bo', 1: 'go', INDETERMINATE: 'ro'}
DEFAULT_COLORS = {0: 'ro', 1: 'go', INDETERMINATE: 'w.'}

#---------------------------------------------
# Methods
#---------------------------------------------
def plot_sequences(haps, title=None, xlabel=None, ylabel=None, x=None, fig=None):
    '''Generate a plot of a list of haplotype sequences/differences d.'''
    xlabel = xlabel if xlabel else 'SNP position' 
    num_sequences = haps.shape[0]
    if fig:
        P.figure(fig)
    P.clf()
    for m in xrange(0, num_sequences):
        P.subplot(num_sequences, 1, m + 1)
        P.ylim([-1.1, 1.1])
        if x is None:
            P.plot(haps[m], 'o')
        else:
            P.plot(x, haps[m], 'o')
        if ylabel is not None:
            P.ylabel(ylabel[m])
    if title is not None:
        P.subplot(num_sequences, 1, 1)
        P.title(title)
    P.subplot(num_sequences, 1, num_sequences)
    P.xlabel(xlabel)

def _xaxis(h, xaxis, snps):
    if type == 'snp':
        return None
    elif xaxis == 'bp':
        return h.snp['base_pair'] / constants.MEGA_BASE_PAIR
    elif xaxis == 'dist_cm':
        return h.snp['dist_cm']

def plot_all_diffs(h, id1, id2, hap1_type=None, hap2_type=None, snps=None, xaxis='snp', fig=None):
    '''Generate a plot of all pairs of hap differences between sample indices id1, id2 in the
    Haplotype object h. xaxis = ''snp'|''bp''.'''
    d = diff.all_diffs(h.data, id1, id2, hap1_type=hap1_type, hap2_type=hap2_type)
    if snps is not None:
        d = d[:, snps]
    else:
        snps = h.snp_range
    plot_sequences(d, title='Hap Differences Between Family Samples %d,%d' % (id1, id2,),
                   x=_xaxis(h, xaxis, snps),
                   # x=haplotype.snp['name'],
                   xlabel='SNP position' if xaxis == 'snp' else 'Mbp',
                   ylabel=['(%d,%d)' % (j, k,) 
                           for j in ((0, 1) if hap1_type is None else (hap1_type,))
                           for k in ((0, 1) if hap2_type is None else (hap2_type,))],
                   fig=fig)
    return d

def plot_ibs_state(g, id1, id2, snps=None, xaxis='snp', fig=None):
    '''Generate a plot of IBS equality between two pairs of genotypes.'''
    gd = g.data
    if snps is not None:
        gd = gd[snps]
    else:
        snps = g.snp_range
    d = diff.ibs_state(gd, id1, id2)
    plot_sequences(np.array([d]), title='IBS Difference Between Family Samples %d,%d' % (id1, id2,),
                   x=_xaxis(g, xaxis, snps),
                      # x=haplotype.snp['name'],
                      xlabel='SNP position' if xaxis == 'snp' else 'Mbp',
                      ylabel=['IBS Difference'],
                   fig=fig)
    P.ylim([-0.3, 2.3])
    return d

def plot_comparison(template, haps, snps=None,
                    title=None, xlabel=None, x=None, y=None, template_y='TEMP',
                    colors=DEFAULT_COLORS):
    '''Generate a recombination coloring plot of the difference between template
    and each element in the list haps..'''
    # Read or re-generate differences data
    num_snps = len(template)
    if np.size(haps.shape) == 1:
        haps = np.transpose(np.array([haps]))
    num_haps = haps.shape[1]
    # d = (haps != np.transpose(np.tile(template, (num_haps,1)))).astype('int8')
    d = np.transpose(np.array([diff.hap_diff(haps[:, j], template).astype('int8') 
                               for j in xrange(0, num_haps)]))
        
    # Restrict view to snps, if specified
    x = x if x is not None else np.arange(0, num_snps)
    if snps is not None:
        (x, d) = (x[snps], d[snps, :])
        
    # Generate haplotype plots
    P.clf()
    P.hold(True)
    hap_ticks = range(0, num_haps)
    for j in hap_ticks:
        difference = d[:, j]
        for (value, color) in colors.iteritems():
            hap = np.where(difference == value)
            if hap:
                P.plot(x[hap[0]], np.zeros((np.size(hap),)) + j, color,
                         markeredgecolor='gray' if value == INDETERMINATE else color[0])
            
    if title is not None:
        P.title(title)
    P.xlim([min(x) - 1, max(x) + 1])
    P.ylim([-0.5, num_haps - 0.5])
    P.xlabel(xlabel if xlabel else 'SNP #')
    P.ylabel('Sample')
    if y is not None:
        P.yticks(hap_ticks, [template_y] + list(y))
    return d

def plot_family_comparison(problem, family, parent_type, template=None, title=None, ylabel=None,
                           colors=DEFAULT_COLORS, snps='het', x=None, y=None, children=None,
                           xaxis=None, yaxis=None):
    '''Generate a recombination coloring plot of the difference between a parent
    and corresponding children haplotype in a nuclear family. Only parent het SNPs are included.
    If template = None, compare against parent; otherwise compare against the template child
    whose ID is template. If a snps array is specified, the plot is restricted to those snps;
    if snp = 'het', all het snps in the parent are used; else, all snps are used.
    
    children: Display all genotyped children (default), or selected children (if specified)'''
    g, h = problem.data
    parent = family.parents[parent_type]
    template_allele = parent_type if template else PATERNAL
    template = template if template is not None else parent
    poo = problem.haplotype.poo_phase
    h_template = h[:, template, 1 - template_allele if poo[template] < 0 else template_allele]  # Template haplotype to use
    # Display all genotyped children, or selected children 
    children = children if children else [child for child in family.children_list if problem.is_genotyped(child)]
    # print parent
    snps = gt.where_heterozygous(g, sample=parent) if snps == 'het' else (snps if snps is not None else problem.snp_range)
    template_y = 'TEMP: %d' % (template,)
    y = y if y is not None else children
    if yaxis == 'id':
        template_y = 'TEMP: %d' % (problem.pedigree.sample_id[template],)
        y = [problem.pedigree.sample_id[child] for child in y]
    ht = h[:, children, parent_type]
    ho = h[:, children, 1 - parent_type]
    flipped = (poo[children] < 0)
    ht[:, flipped] = ho[:, flipped]
    return plot_comparison(template=h_template, template_y=template_y,
                           haps=np.concatenate((h_template[:, np.newaxis], ht), axis=1),
                           snps=snps,
                           xlabel=('SNP position' if xaxis == 'snp' else 'Mbp') if xaxis else None,
                           x=_xaxis(problem.haplotype, xaxis, snps),
                           # x=x if x else ((None if xaxis == 'snp' else problem.haplotype.snp['base_pair'] / constants.MEGA_BASE_PAIR) if xaxis else None),
                           y=y,
                           title='Nuclear Family (%d,%d) Children Comparison: %s Haplotype' % \
                           (family.father, family.mother, ALLELE_LABELS[parent_type].capitalize()),
                           colors=colors), snps,

def plot_all_family_comparisons(problem, family, template=None, title=None, ylabel=None,
                                colors=DEFAULT_COLORS, snps='het', x=None, y=None, xaxis=None, yaxis=None,
                                out_file=None, suffix=ALLELE_LABELS, fig_num=1):
    '''Plot family haplotype comparisons for both parents. Each parent is depicted in a separate
    figure. Save to out_file-prefixed files if specified.'''
    if out_file:
        (file_name, file_ext) = os.path.basename(out_file).split('.')
        dir_name = os.path.dirname(out_file)  
        file_name = (dir_name if dir_name else '.') + '/' + file_name
    d = [None, None]
    for parent_type in ALLELES:
        P.figure(fig_num + parent_type)
        P.clf()
        d[parent_type] = plot_family_comparison(problem, family, parent_type,
                                                template=template, title=title, ylabel=ylabel, colors=colors,
                                                snps=snps, x=x, y=y, xaxis=xaxis, yaxis=yaxis)
        if out_file:
            P.savefig('%s_%s.%s' % (file_name, suffix[parent_type], file_ext))
    return d

def plot_all_family_singlefig(problem, family, children=None,
                              template=None, title=None, ylabel=None,
                              colors=DEFAULT_COLORS, snps='het', x=None, y=None, xaxis=None,
                              out_file=None, suffix=ALLELE_LABELS, fig_num=1):
    '''Plot haplotype comparisons for an entire family or selected children 'children' within a family.
    A single figure is produced with all haplotypes.'''
    raise ValueError('To be implemented')

#---------------------------------------------
# Private Methods
#---------------------------------------------
def plot_hap_corr_matrix(h, members, hap_type=None):
    '''Plot a sib correlation matrix.'''
    # Parametrize all sib pairs
    r = diff.hap_corr_matrix(h, members, hap_type=hap_type)
    P.clf()
    pcolor(r)
    if hap_type:
        P.xlabel('hap %d' % (hap_type[1]))
        P.ylabel('hap %d' % (hap_type[0])) 
    P.title('% Equal Entries Between Haplotypes') 
    P.colorbar()
    return r

def plot_fill_fraction(problem, color='b', zoom=None, ticks=None, label=None):
    '''Plot the filled fraction in each sample of a phased data set.'''
    data = problem.fill_fraction()
    # Sort by ascending fill %
    data = data[np.argsort(data[:, 1])]
    P.plot(range(1, data.shape[0] + 1), data[:, 1], color + '.-', label=label)
    min_y = data[0, 1] * 0.95
    if zoom:
        max_x = np.where(data[:, 1] > zoom)[0][0]
        P.xlim([0, max_x + 1])
    if ticks:
        yticks = linspace(min_y, 1.0, ticks)
        P.yticks(yticks, ['%.3f' % (t,) for t in yticks])
    P.xlabel('Sample')
    P.ylabel('Filled Haplotype %')
    P.grid(True)
    return data

def plot_parametric_experiment(results, print_stats=True):
    '''Plot the results of a parametric validation experiment.'''
    fraction = 100.0 * results['deleted_fraction']
    full_call, = P.plot(fraction, 100.0 * results['full_call_fraction'], 'b-')
    partial_call, = P.plot(fraction, 100.0 * results['partial_call_fraction'], 'b--')
    full_error , = P.plot(fraction, 100.0 * results['full_error_fraction'], 'r-')
    partial_error , = P.plot(fraction, 100.0 * results['partial_error_fraction'], 'r--')
    P.xlabel('% Deleted Genotypes')
    P.ylabel('Rate')
    P.title('Phasing Validation Results')
    P.legend([full_call, partial_call, full_error, partial_error],
               ['Call Rate (Full)', 'Call Rate (Partial)',
                'Error Rate (Full)', 'Error Rate (Partial)'],
               prop={'size': 9})
    P.ylim([-0.01, 100.0])
    if print_stats:
        __print_stat_line('Full    Call Rate' , results['full_call_fraction'])
        __print_stat_line('Partial Call Rate' , results['partial_call_fraction'])
        __print_stat_line('Full    Error Rate', results['full_error_fraction'])
        __print_stat_line('Partial Error Rate', results['partial_error_fraction'])

def plot_experiment_stats(e):
    sample_data = np.where(e.num_test_genotypes(SAMPLE) > 0)[0]
    c_sample = (100.0 * e.called(SAMPLE)[sample_data]) / e.num_test_genotypes(SAMPLE)[sample_data] + 1e-15
    fill = 100.*e.fill[sample_data]

    snp_data = np.where(e.num_test_genotypes(SNP) > 0)[0]
    c_snp = (100.0 * e.called(SNP)[snp_data]) / e.num_test_genotypes(SNP)[snp_data]
    
    # Call % vs. fill %
    P.figure(1);
    P.clf();
    P.plot(fill, c_sample, 'o')
    P.xlabel('Fill %')
    P.ylabel('Call %')
    P.title('Validation Breakdown by Sample, %.2f%% Deleted. r = %.2f' % 
              (100.0 * e.fraction, np.corrcoef(fill + SMALL_FLOAT, c_sample + SMALL_FLOAT)[0, 1],))

    # Call % vs. SNP
    P.figure(2);
    P.clf();
    P.plot(snp_data, c_snp, 'o')
    P.xlabel('SNP #')
    P.ylabel('Call %')
    P.title('Validation Breakdown by SNP, %.2f%% Deleted' % (100.0 * e.fraction,))
    
    return (np.array([snp_data, c_snp]).transpose(),
            np.array([sample_data, c_sample, fill]).transpose())

def pcolor(*args, **kwargs):
    '''Plot a matrix and reverse the y-axis to make it look like a "normal" matrix.'''
    P.pcolor(*args, **kwargs)
    ax = P.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    
#---------------------------------------------
# Manhattan Plots
#---------------------------------------------
def manhattan(subplot, p_values, p=0.05, threshold='bonferroni', title=None, colors=('r', 'g', 'b'),
              min_p_value=1e-16):
    '''Generate a Manahattan plot from a list of 22 p-value data sets. Entry data[i] corresponds
    to chromosome i-1, and should be tuple (snp_bp_coordinate, p_value).'''         
    # Prepare plot
    subplot.set_yscale('log')
    if title:
        P.title(title)
    P.xlabel('Chromosome')
    P.ylabel(r'$p^{-1}$')
    P.hold(True)

    offset = genome_bp_offset()[:NUM_CHROMOSOMES]
    (xmin, xmax) = (np.inf, -np.inf)
    chr_label = []
    for (index, (snp_chr, p_chr)) in enumerate(p_values):
        x = snp_chr + offset[index]
        color = colors[np.mod(index, len(colors))]
        P.scatter(x, 1.0 / np.maximum(min_p_value, p_chr), c=color, marker='o', edgecolor=color)
        xmin = min(xmin, np.min(x))
        xmax = max(xmax, np.max(x))
        chr_label.append('%s' % (index + 1,))
    P.xlim((xmin, xmax))
    # Set the locations and labels of the x-label ticks
    P.xticks(ndimage.convolve(offset, [0.5, 0.5]), chr_label)    
    
    # Calculate significance threshold    
    if threshold == 'bonferroni':
        # Bonferroni correction: divide by the number of SNPs we are considering
        num_snps = sum(len(chr_data[0]) for chr_data in p_values)
        threshold = p / num_snps
    elif isreal(threshold):
        # Custom threshold
        pass
    elif not threshold:
        raise ValueError('Unsupported threshold %s' % (threshold,))
    if threshold is not None:
        # Draw threshold
        P.axhline(y=1.0 / threshold, color='red')

# Deprecated
def prepare_figure_vs_snp(info, title=None, xaxis='snp', snp_index=None):
    '''Prepare a template figure vs. SNP # or location in base pairs for a ProblemInfo object.
    Return a plot data struct.'''
    snp_index = snp_index if snp_index is None else info.snp_range 
    data = util.Struct()
    if xaxis == 'snp':
        data.x = info.snp_range[snp_index]
        data.xlabel = 'SNP #'
    else:
        data.x = info.snp['base_pair'][snp_index] / im.constants.MEGA_BASE_PAIR
        data.xlabel = 'SNP Position [Mbp]'
    data.xmin, data.xmax = np.min(data.x), np.max(data.x)
    if title: P.title(title)
    P.xlabel(data.xlabel)
    P.xlim((data.xmin, data.xmax))
    return data

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __print_stat_line(title, result_field):
    print '%-18s: %7.2f%% +- %.2f%%' % (title, 100 * np.mean(result_field), 100 * np.std(result_field))
    
def genome_bp_offset():
    '''Return the offset required to convert a SNP base pair coordinate to a global genomic
    coordinate. In the returned array, entry i corresponds to chromosome i-1.'''
    return np.concatenate(([0], np.cumsum(db_gene.snp.file_dao.DEFAULT_FILE_DAOS.chrom_dao.total_bp())))
    
def genome_cm_offset():
    '''Return the offset required to convert a SNP base pair coordinate to a global genomic
    coordinate. In the returned array, entry i corresponds to chromosome i-1.'''
    return np.concatenate(([0], np.cumsum(db_gene.snp.file_dao.DEFAULT_FILE_DAOS.chrom_dao.total_cm())))
