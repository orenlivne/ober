'''
============================================================
Generate statistics and plots of IMPUTE2 SNP concordance
with our imputation.

Created on September 17, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib.pyplot as P, numpy as np, util, impute as im, itertools as it, db_gene

class SnpClass(util.Struct):
    def __init__(self, title, (bp, called_in_both, concordance, het_concordance, call_rate, call_rate_pedigree, impute2_info, maf), i):
        '''Filter data to a selected set of SNPs specified by the logical index array i.'''
        super(util.Struct, self).__init__()
        self.title, self.bp, self.called_in_both, self.impute2_info, self.concordance, \
        self.het_concordance, self.maf, self.call_rate, self.call_pedigree = \
        title, bp[i], called_in_both[i], impute2_info[i], concordance[i], het_concordance[i], maf[i], \
        call_rate[i], call_rate_pedigree[i]
    
    def concordance_of_genotype_set(self, genotype_set):
        if genotype_set == 'all': return self.concordance
        elif genotype_set == 'het': return self.het_concordance

    def mean_concordances(self):
        return np.mean(self.concordance), np.mean(self.het_concordance)

def impute2_concordance(in_file):
    '''Calculate impute2 concordance for a single window whose concordance statistics data file is in_file.'''
    # Parameters
    min_called_in_both = 0.4
    maf_threshold = 0.05

    # Columns:
    # 0=x
    # 1=vid
    # 2=bp
    # 3=a1
    # 4=a2
    # 5=#discordant
    # 6=call rate
    # 7=concordance
    # 8=#discordant hets
    # 9=#called hets
    # 10=het_concordance
    # 11=call_rate_impute2
    # 12=call_rate_pedigree
    # 13=impute2_info
    # 14=maf
    usecols = [2, 6, 7, 10, 11, 12, 13, 14]  # Columns #s to use from stats files
    called_in_both_col, maf_col = 1, 7
    
    # TODO: improve maintainability of col #s with a record array with a dtype defining cols+their names
      
    #-----------------------------
    # Load data - phased run 
    #-----------------------------
    data = np.loadtxt(in_file, usecols=usecols) 
    data_cols = tuple(data[:, i] for i in xrange(data.shape[1]))
    # Create SNP classes
    maf, called_in_both = data[:, maf_col], data[:, called_in_both_col]
    all_snps = SnpClass('all', data_cols, (maf > 0) & (called_in_both > min_called_in_both))
    common_snps = SnpClass('common', data_cols, (maf > maf_threshold) & (called_in_both > min_called_in_both))
    # rare_snps = SnpClass('rare', data_cols, (maf <= maf_threshold) & (called_in_both > min_called_in_both))
    return all_snps, common_snps
    
def plot_impute2_concordance((all_snps, common_snps), save_dir=None, plot=False, min_info_to_plot=0.9):
    '''Generate plot of impute2 concordance for a single window from a Struct holding
    statistics on all snps, all_snps.'''
    util.mkdir_if_not_exists(save_dir)
    # Useful variables
    lim_threshold = [0., 1.]
    n = 40
    maf_n = 40
    # info_bins = [0, 0.7, 0.8, 0.85, 0.9, 1]
    info_bins = [0, 0.9, 1]

    threshold = np.linspace(lim_threshold[0], lim_threshold[1], n + 1)
    maf_bins = np.linspace(0, 0.5, maf_n + 1)
    k = 0  # Figure counter

    if save_dir: util.mkdir_if_not_exists(save_dir)
    for snp_class in (all_snps,):  # (all_snps, common_snps, rare_snps): 
#         k += 1
#         P.figure(k)
#         P.clf()
#         P.hold(True)
#         ax = P.subplot(111)
#         ax.scatter(snp_class.concordance, snp_class.impute2_info, color='b', lw=0, label='Info')
#         P.draw()    
#         P.xlabel('Pedigree Imputation-IMPUTE2 Concordance: %s SNPs' % (snp_class.title.capitalize(),))
#         P.ylabel('IMPUTE2 Info')
#         r = np.corrcoef(snp_class.concordance, snp_class.impute2_info)[0, 1]
#         P.title('IMPUTE2 Info vs. Concordance Rate: %s SNPs, r = %.3f' % (snp_class.title.capitalize(), r))
#         P.axvline(x=0.99, linestyle='dotted', linewidth=1, color='k')
#         # P.legend(loc='lower left', prop={'size': 10})
#         P.show()
#         if save_dir: P.savefig(save_dir + '/info-scatter-%s.png' % (snp_class.title,))

        k += 1
        P.figure(k)
        P.clf()
        P.hold(True)
        P.plot(threshold, [np.mean(snp_class.concordance[snp_class.impute2_info >= t]) for t in threshold], 'bo-', label='All Genotypes')
        P.plot(threshold, [np.mean(snp_class.het_concordance[snp_class.impute2_info >= t]) for t in threshold], 'ro-', label='Hets')
        P.legend(loc='lower right', prop={'size': 10})
        P.xlabel('IMPUTE2 Info Threshold')
        P.ylabel('Mean Concordance')
        P.title('Concordance for different IMPUTE2 Threshold: %s SNPs' % (snp_class.title.capitalize(),))
        P.show()
        if save_dir: P.savefig(save_dir + '/mean-concordance-vs-impute2-threshold-%s.png' % (snp_class.title,))

#     # Quality measures along chromosome
#     k += 1
#     P.figure(k)
#     P.clf()
#     P.hold(True)
#     #P.plot(all_snps.bp, all_snps.impute2_info, color='b', label='Info')
#     P.plot(all_snps.bp, all_snps.concordance, color='r', label='Concordance, phased')
#     # if do_gen: P.plot(all_snps.bp, all_snps_gen.concordance, color='g', label='Concordance, unhased')
#     P.legend(loc='lower left', prop={'size': 10})
#     P.xlabel('Base-Pair Position')
#     P.ylabel('$r^2$')
#     P.title('IMPUTE2 Quality Measure Along the Chromosome')
#     P.show()
#     if save_dir: P.savefig(save_dir + '/concordance-vs-bp.png')

    # Plot the mean concordance vs. MAF for different info cut-offs, all genotypes.
    # Contrast with the mean call rate for each cut-off to show that their accuracy-calling
    # trade-off is worse than ours for rare alleles.
    run, data_set = 'phased', all_snps 
    # for run, data_set in (('phased', all_snps),):  # (('phased', all_snps), ('unphased', all_snps_gen)):
    if False: 
        for genotype_set in ('all', 'het'):
            all_concordance = data_set.concordance_of_genotype_set(genotype_set)
            k += 1
            P.figure(k)
            P.clf()
            P.hold(True)
            P.grid(True)
            info_category = np.digitize(data_set.impute2_info, info_bins)
            d = list(it.izip(xrange(1, len(info_bins)), im.plot.colors.get_colors()))
            d = filter(lambda x: info_bins[x[0] - 1] >= min_info_to_plot, it.izip(xrange(1, len(info_bins)), it.repeat('r')))
            for bin_index, color in d:
                high_info = info_category == bin_index
                maf_bin = np.digitize(data_set.maf[high_info], maf_bins)
                concordance = all_concordance[high_info]
    #                call_rate = np.mean(data_set.call_rate[data_set.impute2_info >= info_bins[bin_index - 1]])
                mean_concordance_maf = np.array([(1.*np.mean(concordance[maf_bin == i])) for i in xrange(len(maf_bins))])
                has_data = np.where(np.isfinite(mean_concordance_maf))[0]
    #                 P.plot(maf_bins[has_data], mean_concordance_maf[has_data], color=color, linestyle='solid', marker='o',
    #                        label='info [%.2f,%.2f], call rate %.2f' % \
    #                        (info_bins[bin_index - 1], info_bins[bin_index], call_rate))            
                P.plot(maf_bins[has_data], mean_concordance_maf[has_data], color=color, linestyle='solid', linewidth=2,  # marker='o',
                       label='info [%.2f,%.2f]' % (info_bins[bin_index - 1], info_bins[bin_index]))            
                if len(d) > 1: P.legend(loc='lower right', prop={'size': 10})
                P.xlabel('MAF')
                P.ylabel('Mean Concordance')
                P.title('IMPUTE2 Concordance vs. MAF: %s Run, %s Genotypes (%d SNPs)' % (run.capitalize(), genotype_set.capitalize(), len(all_concordance)))
                P.ylim([0, 1.001])
                P.show()
                if save_dir: P.savefig(save_dir + '/concordance-vs-maf-%s-%s.png' % (run, genotype_set))
                P.ylim([0.95, 1.001])
                if save_dir: P.savefig(save_dir + '/concordance-vs-maf-%s-%s-zoom.png' % (run, genotype_set))

    # Plot the mean concordance vs. MAF for the highest-quality inf category. Compare concordance
    # for all genotypes vs. het genotypes.
    k += 1
    P.figure(k)
    P.clf()
    P.hold(True)
    P.grid(True)
    run, data_set = 'phased', all_snps 
    info_category = np.digitize(data_set.impute2_info, info_bins)
    bin_index = len(info_bins) - 1
    high_info = info_category == bin_index
    maf_bin = np.digitize(data_set.maf[high_info], maf_bins)
    for genotype_set, color in it.izip(('all', 'het'), im.plot.colors.get_colors()):
        all_concordance = data_set.concordance_of_genotype_set(genotype_set)
        concordance = all_concordance[high_info]
        mean_concordance_maf = np.array([(1.*np.mean(concordance[maf_bin == i])) for i in xrange(len(maf_bins))])
        has_data = np.where(np.isfinite(mean_concordance_maf))[0]
        P.plot(maf_bins[has_data], mean_concordance_maf[has_data], color=color, linestyle='solid', linewidth=2,  # marker='o',
                label='%s genotypes' % (genotype_set.capitalize(),))            
    P.legend(loc='lower right', prop={'size': 20})
    P.xlabel('MAF', fontsize=20)
    P.ylabel('Pedigree Imp.-IMPUTE2 Concordance', fontsize=20)
    # P.title('IMPUTE2 Concordance vs. MAF: %s Run, %s Genotypes (%d SNPs)' % (run.capitalize(), genotype_set.capitalize(), len(all_concordance)))
    P.ylim([0, 1.001])
    P.show()
    if save_dir: P.savefig(save_dir + '/concordance-comparison-vs-maf-%s-%s.png' % (run, genotype_set))
    # P.xlim([0.02, 0.5])
    P.ylim([0.97, 1.001])
    if save_dir: P.savefig(save_dir + '/concordance-comparison-vs-maf-%s-%s-zoom.png' % (run, genotype_set))

#    # IMPUTE2 Call Rate
#    k += 1
#    P.figure(k)
#    P.clf()
#    P.grid(True)
#    P.xlim([threshold[0], threshold[-1]])
#    P.bar(threshold - 0.5 * (threshold[-1] - threshold[0]) / n, [(100.*len(all_snps.concordance[all_snps.impute2_info >= t])) / len(all_snps.concordance) for t in threshold], width=(threshold[-1] - threshold[0]) / n)
#    P.xlabel('Minimum Info')
#    P.ylabel('% Genotypes')
#    P.title('IMPUTE2 Call Rate vs. Info Threshold')
#    P.show()
#    if save_dir: P.savefig(save_dir + '/call-rate-vs-impute2_info.png')
#
    k += 1
    P.figure(k)
    P.clf()
    P.hold(True)
    P.grid(True)
    a = np.histogram(all_snps.call_rate, bins=100)
    P.plot(a[1][1:], np.flipud(np.cumsum(a[0][::-1])) / float(sum(a[0])), 'b-', label='All Variants')
    a = np.histogram(common_snps.call_rate, bins=100)
    P.plot(a[1][1:], np.flipud(np.cumsum(a[0][::-1])) / float(sum(a[0])), 'r-', label='Common Variants')
    P.xlabel('Call Rate Threshold (t)')
    P.ylabel('% SNPs with Call Rate >= t')
    P.xlim([0.5, 1])
    P.title('IMPUTE2 Call Rate (Phased)')
    P.legend(loc='lower left', prop={'size': 10})
    P.show()
    if save_dir: P.savefig(save_dir + '/call-rate.png')

def impute2_stats(snp_class):
    a, b = snp_class.mean_concordances()
    return a, b, np.mean(snp_class.bp)
     
def impute2_concordance_of_windows(output_dir, windows, instances_per_node=24):
    '''Return statistics of the windows specified in the tuple array windows.
    Assuming instances_per_node per node slicing.'''
    return [(s,) + impute2_stats(impute2_concordance('%s/node-%04d/run_impute2-%04d.stats.haps' % (output_dir, s / instances_per_node, s))) for s in windows]

def plot_impute2_concordance_chrom(output_dir, chrom, nodes, instances_per_node=24):
    '''Return of all windows along an entire chromosome. Partition is assumed to be into nodes nodes and
    instances_per_node instances in each node.'''
    stats = np.array([(s,) + impute2_stats(impute2_concordance('%s/node-%04d/run_impute2-%04d.stats.haps' % \
                                                               (output_dir, s / instances_per_node, s)))
                                                               for s in xrange(nodes * instances_per_node)],
                     dtype=[
                            ('window', 'i4'),
                            ('concordance', 'f4'),
                            ('het_concordance', 'f4'),
                            ('bp_center', 'i12')
                            ])    
    
    P.figure(1)
    P.clf()
    P.hold(True)
    mbp = stats['bp_center'] / im.constants.MEGA_BASE_PAIR
    P.plot(mbp, stats['concordance'], 'ro-', label='Concordance (all)')
    P.plot(mbp, stats['het_concordance'], 'go-', label='Concordance (hets)')
    P.xlabel('Physical Distance [Mbp]')
    P.ylabel('Quality Measure')
    P.title('IMPUTE2 Performance: Chromosome %d' % (chrom,))
    # P.legend(loc='upper right', prop={'size': 10})
    P.legend(bbox_to_anchor=(0.85, 1), loc=2, borderaxespad=0., prop={'size': 10})
    return stats

def impute2_concordance_chrom_vs_window_size(output_dir_prefix, chrom, runs):
    '''Return the data of all windows along an entire chromosome, for each of several window sizes specified
    via the number of nodes. ''runs'' is an array of tuples (node#s, instances_per_nodes) of each run.'''
    stats = [None] * len(runs)
    for k, (nodes, instances_per_node) in enumerate(runs):
        windows = nodes * instances_per_node
        stats[k] = (nodes, windows, np.array([(s,) + impute2_stats(impute2_concordance('%s.windows_%d/chr%d/run_impute2/node-%04d/run_impute2-%04d.stats.haps' % \
                                                                                       (output_dir_prefix, windows, chrom, s / instances_per_node, s)))
                                          for s in xrange(windows)],
                                          dtype=[
                                                 ('window', 'i4'),
                                                 ('concordance', 'f4'),
                                                 ('het_concordance', 'f4'),
                                                 ('bp_center', 'i12')
                                                 ]))
    result = util.Struct()
    result.chrom = chrom
    result.instances_per_node = instances_per_node
    result.stats = stats
    return result

def plot_impute2_concordance_chrom_vs_window_size(result):
    '''Create plots comparing IMPUTE2 performance for different window sizes.'''
    # Plot concordance vs. window bp location for all genotypes. One plot per window size.
    chrom_total = db_gene.snp.file_dao.DEFAULT_FILE_DAOS.chrom_dao.total_bp_typed()[result.chrom - 1] / im.constants.MEGA_BASE_PAIR
    P.figure(1)
    P.clf()
    P.hold(True)
    for ((_, windows, s), color) in it.izip(result.stats, im.plot.colors.get_colors()):
        w = chrom_total / windows
        mbp = s['bp_center'] / im.constants.MEGA_BASE_PAIR
        P.plot(mbp, s['concordance'], color=color, linestyle='solid', linewidth=2, marker='.', label='w=%.1f Mb (all) mean %.2f' % (w, np.mean(s['concordance'])))
    for ((_, windows, s), color) in it.izip(result.stats, im.plot.colors.get_colors()):
        w = chrom_total / windows
        mbp = s['bp_center'] / im.constants.MEGA_BASE_PAIR
        P.plot(mbp, s['het_concordance'], color=color, linestyle='dashed', linewidth=2, marker='.', label='w=%.1f Mb (het) mean %.2f' % (w, np.mean(s['het_concordance'])))
    P.xlabel('Physical Distance [Mbp]')
    P.ylabel('Concordance')
    P.title('IMPUTE2 Performance on Chromosome %d vs. Avg Window Size w' % (result.chrom,))
    P.legend(loc='lower right', prop={'size': 10})
    # P.legend(bbox_to_anchor=(0.85, 1), loc=2, borderaxespad=0., prop={'size': 10})
