#!/usr/bin/env python
'''
============================================================
Estimate Hutterite kid imputation rates so that Carole can
decide whether to genotype them with a dense or sparse
Illumina chip.

To this end, we select sibs of the target kids, and phase
them at the Illumina SNPs, given the phased imputed haps
of the rest of the Affy people at those SNPs (that are only
75%-called). Then we calculate the coverage of IBD segments
between those kids and the rest of the Affy people to
estimate the sibs' (and therefore the target kids')
imputation call rates.

Created on July 15, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os, numpy as np, matplotlib.pyplot as P, sys
from impute.phasing.phase import NUM_STAGES

def get_child_ids(path):
    '''Load all new Hutt kid IDs from the ID files under the path ''path''.
    Print the subset of kids whose both parents are in the pedigree and genotyped, or just in the pedigree.
    
    Return the list of (kid id, father index, mother index) for all kids whose both parents
    are in the pedigree and genotyped.'''
    m = np.loadtxt(path + '/ids/parents.csv', dtype=int, skiprows=1)
    p = im.hutt_pedigree()
    print 'Untyped kids in Michelle''s study', m.shape[0]
    kids = [x[0] for x in m if p.node_of.has_key(x[1]) and p.node_of.has_key(x[2]) and p.is_genotyped(p.node_of[x[1]]) and p.is_genotyped(p.node_of[x[2]])]
    print 'Kids with genotyped parents', len(kids)
    ids = np.array([(x[0], p.node_of[x[2]], p.node_of[x[1]]) for x in m if p.node_of.has_key(x[1]) and p.node_of.has_key(x[2]) and p.is_genotyped(p.node_of[x[1]]) and p.is_genotyped(p.node_of[x[2]])])
    untyped_parents = [x[0] for x in m if p.node_of.has_key(x[1]) and p.node_of.has_key(x[2]) and not (p.is_genotyped(p.node_of[x[1]]) and p.is_genotyped(p.node_of[x[2]]))]
    print 'Kids with untyped parents that appear in the pedigree', len(untyped_parents), repr(untyped_parents)
    parents_not_in_pedigree = [x[0] for x in m if not (p.node_of.has_key(x[1]) and p.node_of.has_key(x[2]))]
    print 'Kids with whose parents are not both in the pedigree', len(parents_not_in_pedigree), repr(parents_not_in_pedigree)
    return ids

def convert_kids_problem_to_npz(prefix):
    '''Convert PLINK TPED data set to a Problem npz file.'''
    if not os.path.isfile(prefix + '.npz'):
        p = im.io.read_plink(prefix=prefix,
                                   pedigree=im.itu.HUTT_PED, pedigree_genotyped=im.itu.HUTT_GENOTYPED_PED,
                                   haplotype=None, info=None, verbose=True)
        im.io.write_npz(p, prefix + '.npz')

'''Given a list of couples in the pedigree ped, return an array with a kid from each family, if exists.''' 
kid_sibs = lambda ped, couples: np.array(list(set(f.children_list[0] for f in (f for f in (ped.find_family(x[0], x[1]) for x in couples) if f and f.children))))
    
def phase_sibs(prefix, debug):
    '''Phase the sibs as if they were not genotyped, given the haplotype of the rest of the 1415 Affy people.
    Run phasing stages separately in debug mode; otherwise, run the phasing pipeline in one shot.
    Return the final phased Problem object.'''
    stage = NUM_STAGES
    output = prefix + ('.phased' if stage == NUM_STAGES else '.stage%d' % (stage,))
    if os.path.isfile(output + '.npz'):
        # Load existing file if found
        return im.io.read_npz(output + '.npz')
    else:
        # Run phasing
        print 'Phasing sibs'
        if debug:
            for stage in xrange(1, NUM_STAGES + 1):
                output = prefix + ('.phased' if stage == NUM_STAGES else '.stage%d' % (stage,))
                if not os.path.isfile(output + '.npz'):
                    problem = im.phase.main(pedigree=im.itu.HUTT_PED,
                                            input=prefix + ('.stage%d' % (stage - 1,) if stage > 1 else '') + '.npz',
                                            output=output + '.npz', stage=stage, selected_samples=sibs)
        else:
            stage = NUM_STAGES
            output = prefix + ('.phased' if stage == NUM_STAGES else '.stage%d' % (stage,))
            problem = im.phase.main(pedigree=im.itu.HUTT_PED, input=prefix + '.npz', output=output + '.npz',
                                    selected_samples=sibs, min_output=True)
        return problem

def debug_distant_phasing(sample):
    '''Debug distant-phasing of sample ''sample''.'''
    p = im.io.read_npz('/home/oren/ober/out/kids/cytosnp/chr22/cytosnp.imputed.stage4.npz')
    phaser = im.phase_core.new_phaser_chain([im.phase_distant.distant_phaser(single_sample=sample)])
    h = p.haplotype
    print h.fill_fraction(sample=0)
    phaser.run(p, im.PhaseParam())
    print h.fill_fraction(sample=0)

def plot_sib_phasing(chrom, chip, num_children, fill):
    # fill = problem.fill_fraction(sample=sibs)
    # chrom = problem.info.snp['chrom'][0]
    # num_children = [problem.find_family_by_child(s).num_children for s in sibs]
    # P.clf()
    P.figure()
    P.scatter(num_children, fill[:, 1], color='b')
    
    # Add labels
    ax = P.axes()
    ax.set_title('Hutt New Kid Sibs Phasing %%: %s Chip (%d SNPs), Chromosome %d' % \
                 (chip.capitalize(), problem.num_snps, chrom))
    ax.set_xlabel('#Children in Sib Family')
    ax.set_ylabel('Phasing %')

    # P.show()
    
def plot_phasing_errors(chrom, chip, problem, sibs):
    '''Plot the fraction of differences between the original imputed haplotypes of the sibs
    (which we zeroed out) and the rephased haplotypes.'''
    d = im.diff.hap_diff(problem.g, problem.h)
    a = np.concatenate((np.array([(x, problem.genotype.fill_fraction(x)) for x in xrange(problem.num_samples)]), np.sum(d == 1, axis=0, dtype=float) / problem.num_snps), axis=1)
    P.figure(1)
    P.clf()
    P.scatter(a[sibs, 1], a[sibs, 2])
    P.xlabel('Imputation Call Rate')
    P.ylabel('Re-Phasing Error Rate')
    P.title('%s Chip, Chrom %d: Sib Re-Phasing Quality' % (chip, chrom))

def plot_ibd_coverage_chip(chrom, chip, num_children, c):
    '''Generate a plot of sib % IBD coverage for one chip.'''
    # fill = problem.fill_fraction(sample=sibs)
    # chrom = problem.info.snp['chrom'][0]
    # num_children = [problem.find_family_by_child(s).num_children for s in sibs]
    # P.clf()
    P.figure()
    i = np.argsort(num_children)
    print i
    P.plot(c[i, 1], color='b', label='Projected Affy Rate')
    P.plot(c[i, 2], color='r', label='Projected CGI Rate')
    
    # Add labels
    ax = P.axes()
    ax.set_title('Hutt New Kid Sibs IBD Coverage: %s Chip (%d SNPs), Chromosome %d' % \
                 (chip.capitalize(), problem.num_snps, chrom))
    # ax.set_xlabel('#Children in Sib Family')
    ax.set_ylabel('IBD Coverage %')
    ax.legend(loc='lower right', prop={'size': 10})
    # P.show()

def plot_ibd_coverage(chrom, coverage):
    '''Generate a plot of sib % IBD coverage, comparing all chips.'''
    P.figure()
    P.clf()
    P.hold(True)
    P.grid(True)
    chips = coverage.keys()
    i = np.argsort(coverage[chips[0]][:, 3])
    colors = im.plot.colors.fixed
    for k, (chip, c) in enumerate(coverage.iteritems()):
        name = chip.capitalize()
        color = colors[k % len(colors)]
        P.plot(c[i, 1], color=color, marker='D', label='%s - Affy Rate' % (name,))
        P.plot(c[i, 2], color=color, marker='^', label='%s - Allele CGI Rate' % (name,))
        P.plot(c[i, 3], color=color, marker='o', label='%s - Genotype CGI Rate' % (name,))
    
    # Add labels
    ax = P.axes()
    ax.set_title('Hutterite Kids Projected Imputation Call Rates: Chromosome %d' % (chrom,))
    
    # ax.set_xlabel('#Children in Sib Family')
    ax.set_ylabel('IBD Coverage %')
    P.ylim([0, 105])
    ax.legend(loc='lower right', prop={'size': 10})
    # P.show()

def ibd_segments(problem, sib, affy_samples, segment_file):
    '''Calculate IBD segments between sib and the non-sib affy people A.'''
    if not os.path.isfile(segment_file):
        segments = im.ih.between_samples_segments(problem, [sib], affy_samples, params, num_processes=num_processes) 
        segments.save(open(segment_file, 'wb'))
    else: segments = im.segment.SegmentSet.load(open(segment_file, 'rb'))  # Load existing segments
    return segments

'''Total length of a closed segment set.''' 
total_length = lambda segments: sum(x[1] - x[0] + 1 for x in segments)
'''% chrom_length covered by a closed segment set.''' 
coverage_percentage = lambda segments, chrom_length: 100.*float(total_length(im.segment.union(segments))) / chrom_length

def sibs_ibd_coverage(problem, chrom_length, cgi_samples, affy_samples, ibd_coverage_file):
    '''Calculate IBD segment coverage between sibs and A = rest of affy people, which are phased.'''
    print 'Calculating IBD Coverage, file %s' % (ibd_coverage_file,)
    if not os.path.isfile(ibd_coverage_file):
        with open(ibd_coverage_file, 'wb') as f: 
            for sib in sibs:  # [sibs[0]]: #sibs:
                # print 'Calculating sib IBD segments, sample %d, %d processes' % (sib, num_processes)
                segments = ibd_segments(problem, sib, affy_samples, directory + '/segments_%d.out' % (sib,))
                # Calculate IBD coverage (vs. affy people and vs. CGI people) of each sib
                c = ibd_coverage(sib, segments, chrom_length, cgi_samples)
                sys.stdout.write('%d %f %f %f\n' % (sib, c[0], c[1], c[2]))
                f.write('%d %f %f %f\n' % (sib, c[0], c[1], c[2]))
                f.flush()
    return np.loadtxt(ibd_coverage_file)

def ibd_coverage(sib, segments, chrom_length, cgi_samples):
    '''Return the % of the chromosome base pairs covered by the chromosome. chrom_length = length of
    chromosome [bp]. cgi_samples = set of WGS CGI samples.
    Returns a tuple (a,b,c), where
    a = % of genome covered by all segments (~ affy imputation call rate)
    b = % of genome covered by segments shared with CGI samples (~CGI imputation allele call rate)
    c = % of genome covered by segments shared with CGI samples (~CGI imputation allele call rate)
    '''
    # Calculate estimated imputed allele call rate
    cgi_segments = [x.bp for x in segments if set(y[0] for y in x.samples) & cgi_samples]
    # Calculate estimated imputed genotype call rate
    cgi_segments_paternal_hap = [x.bp for x in segments if (sib, im.constants.PATERNAL) in x.samples and set(y[0] for y in x.samples) & cgi_samples] 
    cgi_segments_maternal_hap = [x.bp for x in segments if (sib, im.constants.MATERNAL) in x.samples and set(y[0] for y in x.samples) & cgi_samples]
    cgi_genotype_segments = im.segment.segment_set_op(im.segment.union(cgi_segments_paternal_hap), im.segment.union(cgi_segments_maternal_hap), 'and') 
    return coverage_percentage((x.bp for x in segments), chrom_length), \
coverage_percentage(cgi_segments, chrom_length), \
coverage_percentage(cgi_genotype_segments, chrom_length)


'''Name of data set corresponding to the chip name ''chip''.'''
chip_data_set = lambda chip: 'hutt' if chip == 'affy' else '%s.imputed' % (chip,)

def get_sib_ids(path):
    ped = im.hutt_pedigree()
    ids = get_child_ids(path) 
    sibs = kid_sibs(ped, ids[:, 1:])
    return sibs

read_chip_problem = lambda path, chrom, chip: im.io.read_npz('%s/%s/chr%d/%s.phased.npz' % (path, chip, chrom, chip_data_set(chip)))

####################################################################################
if __name__ == '__main__':

    path = os.environ['OBER_OUT'] + '/kids'
    out_dir = os.environ['OBER'] + '/doc/kids'
    chips = ['cytosnp', 'omniexpress']
    chrom = 22
    debug = 1
    num_processes = 1  # 4
    params = im.PhaseParam()
    
    # Select some genotype sibs of the kids for call rate estimation
    sibs = get_sib_ids(path)
    print 'sibs', sibs
    cgi_samples = set(im.examples.wgs_sample_index())
    
    # for chip in ['cytosnp', 'omniexpress', 'affy']:
    coverage = {}
    for chip in chips:
        print '*' * 80
        print 'Chip', chip
        print '*' * 80
        data_set = chip_data_set(chip)
    
        # Initialize
        directory = path + '/%s/chr%d' % (chip, chrom)
        prefix = '%s/%s' % (directory, data_set)
        print 'Sib Kids Experiment, directory %s' % (directory,)
        convert_kids_problem_to_npz(prefix)
    
        # Phase sibs
        problem = phase_sibs(prefix, False if chip == 'affy' else debug)
        num_snps = problem.num_snps
        print '#SNPs %d' % (num_snps,)
        num_children = [problem.find_family_by_child(s).num_children for s in sibs]
        fill = problem.fill_fraction(sample=sibs)

        # Plot phasing errors
        plot_phasing_errors(chrom, chip, problem, sibs)
        P.savefig('%s/%s_sib_phasing_quality.png' % (out_dir, chip))
        
        # Calculate IBD segments between sibs and A = rest of affy people, which are phased
        bp = problem.info.snp['base_pair']
        chrom_length = bp[-1] - bp[0] + 1
        affy_samples = np.array(list(set(xrange(problem.pedigree.num_genotyped)) - set(sibs)))
        coverage[chip] = sibs_ibd_coverage(problem, chrom_length, cgi_samples, affy_samples, directory + '/ibd_coverage.out')
        
        # Generate plots
        # plot_sib_phasing(chrom, chip, num_children, fill)
        # P.savefig('%s/%s_phasing_fill.png' % (out_dir, chip))
    
    plot_ibd_coverage(chrom, coverage)
    P.savefig('%s/chips_ibd_coverage.png' % (out_dir,))
