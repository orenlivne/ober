'''
============================================================
Impute SNPs from 98 CGI WGS samples to all Hutt samples
using an identity-by-state segment dictionary (a 
SmartSegmentSet).

Created on February 2, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, os, matplotlib.pyplot as P, util
from impute.phasing.examples import OUT_PHASING
from db_gene.cgi.ids import DEFAULT_ID_INDEX
from impute.data.problem import Problem
from impute.imputation.reader import problem_to_imputation_set
from impute.cgi.extract_genotypes import extract_genotypes
from impute.imputation.ImputationSet import ImputationSet
from impute.data.constants import CHROMOSOMES

#---------------------------------------------
# Methods
#---------------------------------------------
def impute_problem(problem, index_file=DEFAULT_ID_INDEX, samples=None, snps=None, debug=0,
                   input_dir=None, algorithm='index', segment_location=None, ibd_sample_index=None,
                   remove_partial_calls=False, debug_sample=-1, genotypes_as_is=False):
    '''Impute the SNPs that have also been genotyped for all samples, and compare with the
    imputation result. Returns the phasing Problem and the imputed results.
    
    problem = npz file location or Problem instance.
    segment_location = IBD segment repository (segment file or segment index) location.
    ibd_sample_index = optional dictionary of CGI-ID-problem-sample-ID. If None, the mapping between the two IDs is assumed to be the identity.
    
    algorithm options:
    'index' - Imputation V2 using IBD cliques at each SNP (default).
    'interval_tree' - Imputation V1 using interval tree queries on a SmartSegmentSet.'''
    # Read input data
    p = problem if isinstance(problem, Problem) else im.io.read_npz(problem)
    
    # Select SNPs
    if snps is not None: p = p.sub_problem_of_snps(snps)    
    t = problem_to_imputation_set(p, index_file=index_file, genotypes_as_is=genotypes_as_is)
    print p
    print t.genotype
    
    # Run imputation
    segment_location = segment_location if segment_location else ('%s/segments.out' % (input_dir,) if algorithm == 'interval_tree' 
                                                                  else '%s/index_segments' % (os.environ['OBER_OUT'],))
    print 'IBD segment index at', segment_location    
    print 'Imputing, algorithm %s ... ' % (algorithm,)
    if algorithm == 'interval_tree':
        ibd = im.smart_segment_set.SmartSegmentSet.load(p.pedigree.num_genotyped, segment_location)
        print 'IBD segments', ibd.size
        im.imputation.iibd.impute(p.haplotype, ibd, t, samples=samples, debug=debug, ibd_sample_index=ibd_sample_index)
    elif algorithm == 'index':
        ibd = im.index.segment_index.SegmentIndex(segment_location)
        im.imputation.impute_ibd_index.impute(ibd, t, samples=samples, debug=debug, genotype=p.g, ibd_sample_index=ibd_sample_index, remove_partial_calls=remove_partial_calls, debug_sample=debug_sample)
    else:
        raise ValueError('Unsupported imputation algorithm ''%s''' % (algorithm,))
    return p, t

def impute_chrom(chrom=None, snp_step_size=None, snps=None, debug=0, input_dir=None, algorithm='index',
                 segment_location=None, remove_partial_calls=False):
    '''Impute a phased chromosome. Uses default file locations. Imputed SNPs that are uniformly
    spaced (in snp index) with step size snp_step_size.
    Returns Problem, ImputationSet, snp range objects.'''
    if chrom:
        input_dir = input_dir if input_dir else '%s/chr%d' % (OUT_PHASING, chrom)
        problem = im.io.read_npz('%s/hutt.phased.npz' % (input_dir,))
        snps = snps if snps is not None else (np.arange(0, problem.num_snps, snp_step_size) if snp_step_size else None)
    else:
        if not snp_step_size: raise ValueError('Must specify a positive snp_step_size if loading entire genome')
        problem = load_genome_problem(input_dir, snp_step_size)
        snps = None
    p, t = impute_problem(problem, input_dir=input_dir, algorithm=algorithm, segment_location=segment_location,
                          snps=snps, debug=debug, remove_partial_calls=remove_partial_calls)
    return p, t

stats_computer = lambda t, g: im.imputation.istat.StatsComputer(t, g)

def plot_stats(stats, save_prefix=None, fig_num=1, snp_style='continuous',
               filter_length=20):
    '''Plot imputation statistics for a phased chromosome, validating against the original genotypes.'''
    if save_prefix: util.mkdir_if_not_exists(os.path.dirname(save_prefix))
     
    P.figure(fig_num)
    P.clf()
    # P.show()
    stats.plot_vs_snp(snp_style=snp_style, x_axis='cm_cumulative', filter_length=filter_length)
    if save_prefix: P.savefig(save_prefix + '-snp-cm-cumulative.png')  

    fig_num += 1
    P.figure(fig_num)
    P.clf()
    # P.show()
    stats.plot_vs_snp(snp_style=snp_style, x_axis='bp_cumulative', filter_length=filter_length)
    if save_prefix: P.savefig(save_prefix + '-snp-bp-cumulative.png')  

    fig_num += 1
    P.figure(fig_num)
    P.clf()
    # P.show()
    stats.plot_vs_snp(snp_style=snp_style, x_axis='cm_edge_dist', filter_length=filter_length)
    if save_prefix: P.savefig(save_prefix + '-snp-cm-edge-dist.png')  

    fig_num += 1
    P.figure(fig_num)
    P.clf()
    # P.show()
    stats.scatter_snp_concordance(snp_style=snp_style)
    if save_prefix: P.savefig(save_prefix + '-snp-concordance.png')

    fig_num += 1
    P.figure(fig_num)
    P.clf()
    # P.show()
    stats.plot_vs_maf(snp_style=snp_style)
    if save_prefix: P.savefig(save_prefix + '-maf.png')  

    fig_num += 1
    P.figure(fig_num)
    P.clf()
    # P.show()
    stats.plot_vs_sample()
    if save_prefix: P.savefig(save_prefix + '-sample.png')

    if save_prefix: stats.summary(open(save_prefix + '-stats.txt', 'wb'))        

def pipeline_impute_chrom(chrom=None, snp_step_size=None, snps=None, debug=0, input_dir=None, algorithm='index',
                          sample_call_rate_threshold=0.5, plot=False, save_location=None, segment_location=None):
    '''Runs the entire Affymetrix chromosome + stats validation pipeline.
    plot: generate plots if and only if this flag is true.
    save_location: save plots as PNG files under this location. Only relevant if plot=True.'''
    p, t = im.v.iv.impute_chrom(chrom, snp_step_size=snp_step_size, snps=snps, debug=debug,
                                input_dir=input_dir, algorithm=algorithm, segment_location=segment_location)
    computer = stats_computer(t, p.g)
    stats = computer.stats()
    stats.summary()
    if plot: plot_stats(stats, save_prefix=save_location + '/all' if save_location else None, fig_num=1)

    # Remove samples with low call rates and recalculate statistics. Threshold suggested by Gorka
    # on 3-APR-2013 to have a concordance > 98% and look like a round number for a paper publication.
    print ''
    samples = stats.samples_with_call_rates_above(sample_call_rate_threshold)
    stats2 = computer.stats(samples=samples)
    stats2.summary()
    if plot: plot_stats(stats2, save_prefix=save_location + '/selected' if save_location else None, fig_num=4)
    return p, t, computer, stats, stats2

def pipeline_validation_experiment(location_file, true_type, true_location, pedigree, debug=False, remove_partial_calls=False):
    '''Load (the ''true'') genotypes from an external source. Load a list of locations from ''location_file''. Impute them and compare
    with the true genotypes.'''
    g = extract_genotypes(location_file)
    t = ImputationSet(pedigree, g)
    if true_type == 'iplex': true_genotype = im.imputation.reader.iplex_to_genotype(true_location, t)  # os.environ['OBER'] + '/data/impute/rare/to_livne_20121205', t)
    else: raise ValueError('Unsupported true genotype format ''%s''' % (true_type,))
    problem = Problem(pedigree, true_genotype)
    p, t = impute_problem(problem, debug=debug, remove_partial_calls=remove_partial_calls)
    return p, t

def pipeline_validation_hutt_plink(data, out, debug=False, snp_style='continuous', save_location=None,
                                   sample_call_rate_threshold=0.5, filter_length=20):
    '''Run a validation experiment on the PLINK TPED data set ''data'' (=plink data set prefix).
    The data set must match the ordering of the standard Hutterites pedigree, im.itu.HUTT_PED.
    Save results to the output directory ''out''. If ''save_location'' is specified, generate performance
    plots and save them under this directory'''
    imputation_set_file = out + '/imputation_set.npz'
    stats_file = out + '/stats.npz'
    stats_all_file = out + '/stats-all.npz'
    stats_selected_file = out + '/stats-selected.npz'
    if os.path.exists(stats_file) and os.path.exists(stats_all_file) and os.path.exists(stats_selected_file):
        stats = im.imputation.impute_stats.Stats.load(stats_all_file)
        stats_all = im.imputation.impute_stats.Stats.load(stats_all_file)
        stats_selected = im.imputation.impute_stats.Stats.load(stats_selected_file)
    else:
        # Convert TPED to Problem npz if not already done
        if os.path.exists(data + '.npz'): problem = im.io.read_npz(data + '.npz')
        else:
            problem = im.io.read_plink(pedigree=im.itu.HUTT_PED, pedigree_genotyped=im.itu.HUTT_GENOTYPED_PED, 
                                       genotype=data + '.tped', debug=True, frames=None, lam=None)
            im.io.write_npz(problem, data + '.npz')

        # Impute and save results to ImputationSet npz file if not already done
        if os.path.exists(imputation_set_file): p, t = problem, im.v.iv.ImputationSet.load(imputation_set_file)
        else:
            p, t = im.v.iv.impute_problem(problem, debug=debug)
            t.save(imputation_set_file)

        # Compute and save stats in an npz file if not already done
        print 'Computing stats...'
        computer = im.v.iv.stats_computer(t, p.g)
        stats = computer.stats()
        stats.save(stats_file)

        # Remove a sample with very low call rate, likely due to sample contamination since it is not
        # IBS>=1 with its father only in 97% of the genome (should be 99%+ as for other parent-child pairs)
        samples = stats.samples_with_call_rates_above(0.1)
        stats_all = computer.stats(samples=samples)
        stats_all.save(stats_selected_file)

        # Remove samples with low call rates and recalculate statistics. Threshold suggested by Gorka
        # on 3-APR-2013 to have a concordance > 98% and look like a round number for a paper publication.
        samples = stats.samples_with_call_rates_above(sample_call_rate_threshold)
        stats_selected = computer.stats(samples=samples)
        stats_selected.save(stats_selected_file)

    # Plot statistics
    stats_all.summary()
    im.v.iv.plot_stats(stats_all, save_prefix=save_location + '/all' if save_location else None, fig_num=1,
                       snp_style=snp_style, filter_length=filter_length)
    print ''
    stats_selected.summary()
    im.v.iv.plot_stats(stats_selected, save_prefix=save_location + '/selected' if save_location else None, fig_num=100,
                       snp_style=snp_style, filter_length=filter_length)

    return stats, stats_all, stats_selected
    
#---------------------------------------------
# Methods - Merge Entire Geome
#---------------------------------------------
def load_chrom_problem(chrom, input_dir=None, snp_step_size=None):
    '''Load a chromosomal problem. Down-sample SNPs by snp_step_size.'''
    if not input_dir: input_dir = im.examples.OUT_PHASING
    p = im.io.read_npz('%s/chr%d/hutt.phased.npz' % (input_dir, chrom))
    return p.sub_problem_of_snps(np.arange(0, p.num_snps, snp_step_size) if snp_step_size else None)

def load_genome_problem(input_dir=None, snp_step_size=None):
    '''Load the entire genome (merge all chromosomal problems). Down-sample SNPs by snp_step_size.'''
    return sum(im.v.iv.load_chrom_problem(c, snp_step_size=snp_step_size) for c in CHROMOSOMES)
