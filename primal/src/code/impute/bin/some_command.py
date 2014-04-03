#!/usr/bin/env python
'''
============================================================
Scratch space to run and debug commands.

Created on August 16, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, impute as im, impute.batch.batch_util as bu, numpy as np, cProfile, pstats, os, util, pickle, matplotlib.pyplot as P
from impute.phasing.examples import CHR22, OUT_PHASING
from impute.imputation.ImputationSet import ImputationSet
import linecache

#----------------------------------------------------------------------------
# Profiling identity coefficient calculation
#ibd = im.index.segment_index.SegmentIndex(os.environ['OBER'] + '/out/index_segments')

# def idc():
#     t_counter, total_bp = im.poo.idcoef.idcoef(ibd, 22, out='count', do_region=0)
#     
# cProfile.run('idc()', 'idc')
# p = pstats.Stats('idc').strip_dirs()
# p.sort_stats('cumulative').print_stats(50)

#t_counter, total_bp = im.poo.idcoef.idcoef(ibd, 22, algorithm='full', out='count', do_region=0)

#----------------------------------------------------------------------------
# Non-QF PO distortion validation 
# p, t = im.v.iv.impute_chrom(22, snps=[409], debug=1, input_dir=os.environ['OBER_OUT'] + '/phasing/chr22')

#----------------------------------------------------------------------------
# problem = im.io.read_plink(prefix=os.environ['OBER_OUT'] + '/requests/monogenic/work/monogenic.12', pedigree=im.itu.HUTT_PED, haplotype=None, frames=None)
# p, t = im.imputation.qc.qc_problem(problem, debug=2)

#----------------------------------------------------------------------------
# Validation on non-QC affy SNPs
# Not needed should be the same
# p = im.io_pedigree.read('c:/ober/data/hutt/hutt.3chips.noqc.downsampled.tfam', genotyped_id_file='c:/ober/testdata/pedigree/hutterites.tfam')
# p = im.io.read_plink(pedigree=im.itu.HUTT_PED, genotype=os.environ['OBER_DATA'] + '/hutt/hutt.3chips.noqc.downsampled.tped', debug=True, frames=None, lam=None)

#----------------------------------------------------------------------------
# line = linecache.getline('c:/ober/out/impute_cgi/impute2.windows_24/chr22/run_impute2/test.tped', 1)
# print line

# s = im.v.impute2.impute2_validation.impute2_concordance(os.environ['OBER_OUT']+'/impute_cgi/impute2/impute2.stats.clean.downsampled2')

# im.v.impute2.impute2_validation.plot_impute2_concordance_chrom_vs_window_size(os.environ['OBER_OUT'] + '/impute_cgi/impute2', 22, [1, 2, 4])

#----------------------------------------------------------------------------
# problem = im.io.read_plink(prefix=os.environ['OBER_OUT'] + '/kids/reimputation/edge-snp', pedigree=im.itu.HUTT_PED, haplotype=None, frames=None)
# p, t = im.v.iv.impute_problem(problem, debug=2)
# p = im.hutt('hutt.phased.npz')
# im.io_genotype.write('impute2', p.haplotype, sys.stdout, samples=np.arange(5), snps=np.array([1565]), flip_alleles=np.array([1], dtype=bool))

#----------------------------------------------------------------------------
# Sub-problem poo phase saving
# p = im.hutt('hutt.stage3.npz')
# q = p.sub_problem_of_parents(1822, 1825)
# print q.haplotype.poo_phase

#----------------------------------------------------------------------------
# Debugging IBD segments of Hutt kids
# im.examples.hutt_ibd_segments('hutt.phased.npz', 1, None, 164, None, debug=True, min_len=1, len_unit='cm')
# problem = im.io.read_npz('/home/oren/ober/out/kids/cytosnp/chr22/cytosnp.imputed.phased.npz')
# s = im.examples.problem_ibd_segments(problem, 465, 0, 677, 1, debug=True)

#----------------------------------------------------------------------------
# Monogenics validation
# problem = im.io.read_plink(prefix=os.environ['OBER_OUT'] + '/requests/monogenic/work/monogenic.12', pedigree=im.itu.HUTT_PED, haplotype=None, frames=None)
# p, t = im.v.iv.impute_problem(problem, debug=2, remove_partial_calls=True)
#
# im.io.write_plink(im.Problem(genotype=t.imputed, pedigree=im.examples.hutt_pedigree(), haplotype=None, frames=None),
#                  os.environ['OBER_OUT'] + '/requests/monogenic/work/imputed.12', save_frames=False, save_haplotype=False)

# p, t = im.v.iv.pipeline_validation_experiment(open(os.environ['OBER_OUT'] + '/requests/rare/iplex-snp-coords.txt', 'rb'), 'iplex', os.environ['OBER_OUT'] + '/requests/rare/to_livne_20121205', im.examples.hutt_pedigree(), debug=2)

# p, t = im.v.iv.impute_chrom(22, snp_step_size=200, debug=2, input_dir=os.environ['OBER_OUT'] + '/phasing/chr22')
# t.save()
# p = im.io.read_npz(os.environ['OBER_OUT'] + '/imputed/chr22-problem.npz')
# t = ImputationSet.load(os.environ['OBER_OUT'] + '/imputed/chr22-imputed.npz')
# stats = im.v.iv.stats_computer(t, p.g).stats()
# stats.plot_vs_snp(snp_style='discrete')
# P.show()

#----------------------------------------------------------------------------
# def index():
#     im.index.index_segments.main(None, os.environ['OBER_OUT'] + '/phasing/chr12' + '/hutt.phased.info.npz',
#                                  os.environ['OBER_OUT'] + '/phasing/chr12/segments.out',
#                                  '-',
#                                  min_len=0.0, debug=1, algorithm='amg',
#                                  num_processes=1, region_size=10, regions='350')
# 
# cProfile.run('index()', 'index')
# p = pstats.Stats('index').strip_dirs()
# p.sort_stats('cumulative').print_stats(50)

#----------------------------------------------------------------------------
# p, t, i = im.v.iv.impute_chrom(22, snps=[3000], debug=2, input_dir=os.environ['OBER_OUT'] + '/phasing/chr22', segment_file=os.environ['OBER_OUT'] + '/ibd/chr22/segments.3000.out', algorithm='interval_tree')

# im.examples.hutt_ibd_segments('hutt.phased.npz', 81, 0, 987, 1, debug=True)
# im.index.index_segments.main(None, im.examples.CHR22 + '/hutt.phased.info.npz', im.examples.CHR22 + '/segments.out', os.environ['OBER_OUT'] + '/ibd/chr22/index_segments_test', num_processes=1, region_size=1, snp_index=456, min_len=0.0, debug=2, algorithm='amg')
# p, t, i = im.v.iv.impute_chrom(22, snps=[456], debug=2, input_dir= os.environ['OBER_OUT'] + '/phasing/chr22', segment_file=os.environ['OBER_OUT']+'/ibd/chr22/index_segments_test', algorithm='index')

#----------------------------------------------------------------------------
# p, t, i, s, stat = im.v.iv.impute_chrom_pipeline(1, snps=[20400], debug=True, input_dir=OUT_PHASING + '/chr1', save=False)
# p, t, i = im.v.iv.impute_chrom(22, snps=[3060], debug=2, input_dir=OUT_PHASING + '/chr22')

# def imp():
#    p, t, i = im.v.iv.impute_chrom(22, snps=[3060], debug=2, input_dir=OUT_PHASING + '/chr22')
# cProfile.run('imp()', 'imp')
# p = pstats.Stats('imp').strip_dirs()
# p.sort_stats('cumulative').print_stats(50)

#----------------------------------------------------------------------------
# def sim():
#    call_rate_sim(util.Struct(debug=True, num_simulations=50, num_alphas=1, ped_file='/home/oren/ober/testdata/pedigree/hutterites.ped'))
#
# cProfile.run('sim()', 'sim')
# p = pstats.Stats('sim').strip_dirs()
# p.sort_stats('cumulative').print_stats(50)
#----------------------------------------------------------------------------
# im.examples.hutt_ibd_segments('hutt.phased.npz', 1007, 1, 1041, 1, debug=True)
# im.examples.hutt_ibd_segments('hutt.phased.npz', 997, 0, 1041, 1, debug=True)
# im.examples.hutt_ibd_segments('hutt.phased.npz', 997, 0, 1007, 1, debug=True)
# im.examples.hutt_ibd_segments('hutt.phased.npz', 1220, 1, 1007, 1, debug=True)
# im.examples.hutt_ibd_segments('hutt.phased.npz', 0, 0, 188, 0, debug=True)
# im.examples.hutt_ibd_segments('hutt.phased.npz', 11, 0, 281, None, debug=True)
# i = [2720] 
# t = im.examples.impute_chr22(i, debug=True)

#----------------------------------------------------------------------------
# def impute():
#    p = im.hutt('hutt.phased.npz')
#    ibd = im.smart_segment_set.SmartSegmentSet.load(p.pedigree.num_genotyped, os.environ['OBER_OUT'+'/segments.out')
#    t = im.imputation.ImputationSet.from_file(p.pedigree, im.itu.IMPUTE_RARE_SAMPLE) #@UndefinedVariable
#    snps = [np.where(t.snp['chrom'] == 22)[0][0]] # SNP list out of all SNPs in t to impute 
#    im.imputation.iibd.impute(p.haplotype, ibd, t, snp=snps, debug=False)
#
# def impute_test():
#    p = im.io.read_npz(im.itu.IMPUTE_PHASED_SAMPLE)
#    ibd = im.smart_segment_set.SmartSegmentSet.load(im.itu.IMPUTE_IBD_SAMPLE)
#    t = im.imputation.ImputationSet.from_file(p.pedigree, im.itu.IMPUTE_RARE_SAMPLE) #@UndefinedVariable
#    snps = np.where(t.snp['chrom'] == 22)[0] # SNP list out of all SNPs in t to impute 
#    im.imputation.iibd.impute(p.haplotype, ibd, t, snp=snps, debug=False)
#
# cProfile.run('impute()', 'impute')
# p = pstats.Stats('impute').strip_dirs()
# p.sort_stats('cumulative').print_stats(50)

# cProfile.run('impute_test()', 'impute_test')
# p = pstats.Stats('impute_test').strip_dirs()
# p.sort_stats('cumulative').print_stats(50)

###----------------------------------------------------------------------------
# def stage0():
#    im.examples.phase_chr22('hutt.npz', 'hutt.phased.npz', stage=0) #, debug=True)
#
# cProfile.run('stage0()', 'stage0')
# p = pstats.Stats('stage0').strip_dirs()
# p.sort_stats('cumulative').print_stats(50)

#----------------------------------------------------------------------------
# p = im.hutt('hutt.stage2.npz')
# a = p.sub_problem_of_parents(7, 462)
# print a.info.ibd
# a = p.sub_problem_of_parents(79, 324)
# print a.info.ibd

#----------------------------------------------------------------------------
# p = im.io.read_npz(im.itu.CHR18_FAMILY_ISOLATED)
# phaser = im.phase_main.main_phaser(print_times=True)
# phaser.run(p, im.PhaseParam(debug=True))
# f = list(p.families(genotyped=False))[2]
# d = im.plots.plot_all_family_comparisons(p, f, template=0)

#----------------------------------------------------------------------------
# phasing_dir = os.environ['OBER'] + 'phasing.20121130'
# reader = PhasingResultReader(phasing_dir, 'hutt_phased')
# t = ImputationSet(reader.pedigree, os.environ['OBER'] + '/data/impute/rare/rare.npz')
# q = im.imputation.reader.iplex_to_genotype(os.environ['OBER'] + '/data/impute/rare/to_livne_20121205', t)
# print q.data[0, t.sample_index, :]
# print t.genotype.data[0, :, :]
# print q.data[11, t.sample_index, :] - t.genotype.data[11, :, :]
# pass

##----------------------------------------------------------------------------
# def first_nz(d):
#    for i in xrange(len(d)):
#        if d[i] != 0:
#            return i
#
# import time
# a = np.random.random_integers(0, 10, 100)
# t0 = time.time()
# for i in xrange(100):
#    (a != 0).nonzero()[0][:]
# t1 = time.time()
# print (t1-t0)/100

#----------------------------------------------------------------------------
# chrom = 1
# p = np.load('phasing.20121129/hutt_phased_chr%d.stats.npz' % (chrom,))
# ibd = im.ibd.segment.SegmentSet(sum((x.ibd for x in p['info'][0]), []))

#----------------------------------------------------------------------------
# np.set_printoptions(linewidth=200)
#
# p = im.hutt('hutt.stage1.npz')
# print p.h[2955:2965, 370, 1]
# im.phase.family_phaser().run(p, im.PhaseParam(debug=True, snp_test_index=2960))
# print p.h[2955:2965, 370, 1]

# p2 = p.sub_problem(p.pedigree.neighbors(1298, depth=4, genotyped=False))
# print p2.info.genotype_frequency.shape

#----------------------------------------------------------------------------
# p = im.examples.phase_chr22('hutt.npz', 'hutt.stage1.npz', stage=1)

#----------------------------------------------------------------------------
# p = im.hutt('hutt.stage4.npz')
# p2 = p.sub_problem(p.pedigree.neighbors(1298, depth=4, genotyped=False))
# print p2.info.genotype_frequency.shape

##----------------------------------------------------------------------------
# #p = im.io.read_npz(im.itu.FAMILY13_STAGE3)
# p = im.io.read_npz(im.itu.FAMILY4_STAGE3)
# m = im.ig.ibd_germline(p, p.first_family.member_list)
# m.pprint_segments(True)
# d = m.group_by_segment()
# print d

#----------------------------------------------------------------------------
# p = im.hutt('hutt.stage2.npz')
# f = p.find_family(1503, 945, genotyped=False)
# i=2829 #507 #945
# f = p.find_families_by_member(i, genotyped=False).next()
# phaser.handle(util.Struct(problem=p, params=im.PhaseParam(debug=True, single_member=i)))
# p = im.io.read_npz(im.itu.FAMILY945_ONE_PARENT_STAGE2)
# f = p.families(genotyped=False)[0]
# print p.fill_fraction(sample=f.member_list)
# phaser = im.phase.family_child_comparison_phaser()
# phaser.handle(util.Struct(problem=p, params=im.PhaseParam(debug=True)))
# #phaser.handle(util.Struct(problem=p, params=im.PhaseParam(debug=True)))
# print p.fill_fraction(sample=f.member_list)
#----------------------------------------------------------------------------
# #examples.generate_test_sets()
# p = im.hutt('hutt.stage4.npz')
# im.io.write_npz(p.sub_problem_of_family(p.find_families_by_father(1715, genotyped=False)[0], im.itu.FAMILY963_STAGE4))
#
# p = im.hutt('hutt.npz')
# phaser = im.phase_main.main_phaser()
# e = v.Experiment(p, 0.01)
# e.run(phaser)
