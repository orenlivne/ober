#!/usr/bin/env python
'''
============================================================
Examples of running main phasing scripts for use at the
python command line.

Created on July 26, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, impute as im, numpy as np, itertools, db_gene

#---------------------------------------------
# Constants
#---------------------------------------------
PREFIX = os.environ['OBER_OUT']
OUT_PHASING = PREFIX + '/phasing'
CHR22 = OUT_PHASING + '/chr22'

#---------------------------------------------
# Methods
#---------------------------------------------
def test():
    print 'hi'

def phase_chr22(infile, outfile, stage=0, debug=False, **kwargs):
    '''Run phasing on chromosome 22.'''
    return im.phase.main(pedigree=im.itu.HUTT_PED, input=CHR22 + '/' + infile,
                         output=CHR22 + '/' + outfile, stage=stage, debug=debug, **kwargs)
    
'''Impute validation on a set of chroosome 22snp indices.'''
impute_chr22 = lambda snps, samples = None, debug = False: im.v.iv.impute_problem(CHR22 + '/hutt.phased.npz', samples=samples, snps=snps, debug=debug)

def npz_ibd_segments(npz_file, i, ai, j, bj, **kwargs):
    '''IBD segments using a phasing npz file- absolute path.'''
    return problem_ibd_segments(im.io.read_npz(npz_file), i, ai, j, bj, **kwargs)

def hutt_ibd_segments(hutt_file, i, ai, j, bj, **kwargs):
    '''IBD segments using a phasing npz file relative to the chr22 directory.'''
    return problem_ibd_segments(im.hutt(hutt_file), i, ai, j, bj, **kwargs)

def problem_ibd_segments(problem, i, ai, j, bj, **kwargs):
    '''IBD segments using a phasing npz file- absolute path.'''
    segments = im.segment.SegmentSet()
    for a in (im.constants.ALLELES if ai is None else [ai]):
        for b in (im.constants.ALLELES if bj is None else [bj]):
            segments += im.ih.hap_segments_from_problem(problem, (i, a), (j, b), im.PhaseParam(**kwargs))
    return segments

'''Hutterites Pedigree object.'''
hutt_pedigree = lambda _ : im.io_pedigree.read(im.itu.HUTT_PED, genotyped_id_file=im.itu.GENOTYPE_SAMPLE + '.tfam')

'''Load full hutt problem in npz format from chr22 data directory.'''
hutt = lambda npz_file: im.io.read_npz(CHR22 + '/' + npz_file)

def hutt_plink(stage):
    '''Load full hutt problem in plink format from chr22 data directory.'''
    return im.io.read_plink(pedigree=im.itu.HUTT_PED,
                            prefix=CHR22 + '/hutt.' + stage,
                            pedigree_genotyped=CHR22 + '/hutt.tfam')
    
def hutt_pedigree():
    '''Load hutt pedigree.'''
    return im.io_pedigree.read(im.itu.HUTT_PED, genotyped_id_file=CHR22 + '/hutt.tfam')

def generate_test_sets():
    '''Re-generate all chr22 test sets.'''     
    generate_test_sets_phasing()
    generate_test_sets_segments()

def generate_test_sets_phasing():
    '''Re-generate all chr22 phasing test sets.'''     
    p = hutt('hutt.stage1.npz')
    im.io.write_npz(p.sub_problem_of_parents(225, 81), im.itu.FAMILY225_STAGE1)
    im.io.write_plink(p.sub_problem_of_parents(68, 184), im.itu.FAMILY7)
    im.io.write_npz(p.sub_problem_of_parents(79, 324), im.itu.FAMILY_TOO_ZEROED_STAGE1)
    im.io.write_npz(p.sub_problem_of_parents(505, 533), im.itu.FAMILY4_STAGE1)
    
    p = hutt('hutt.stage2.npz')
    im.io.write_npz(p.sub_problem_of_parents(7, 462), im.itu.FAMILY13_STAGE2)
    f = p.find_families_by_parent(im.constants.MATERNAL, 12)[0]
    im.io.write_npz(p.sub_problem_of_parents(f.father, f.mother), im.itu.FAMILY12_STAGE2)
    im.io.write_npz(p.sub_problem_of_parents(225, 81), im.itu.FAMILY225_STAGE2)
    im.io.write_npz(p.sub_problem_of_parents(1503, 945), im.itu.FAMILY945_ONE_PARENT_STAGE2)
    im.io.write_npz(p.sub_problem_of_parents(79, 324), im.itu.FAMILY_TOO_ZEROED_STAGE2)

    p = hutt('hutt.stage3.npz')
    im.io.write_npz(p.sub_problem_of_parents(7, 462), im.itu.FAMILY13_STAGE3)
    im.io.write_plink(p.sub_problem_of_parents(7, 462), im.itu.FAMILY13_STAGE3_PLINK)
    im.io.write_npz(p.sub_problem_of_parents(2003, 1415), im.itu.FAMILY2003_STAGE3)
    im.io.write_npz(p.sub_problem_of_parents(505, 533), im.itu.FAMILY4_STAGE3)
    im.io.write_npz(p.sub_problem_of_parents(1532, 1924), im.itu.SIB_FOUNDERS_STAGE3)
    im.io.write_npz(p.sub_problem_of_parents(1822, 1825), im.itu.SIB_FOUNDERS_1049_STAGE3)

    p = hutt('hutt.stage4.npz')
    im.io.write_npz(p.sub_problem(p.pedigree.neighbors(1298, depth=4, genotyped=False)), im.itu.NBHRS1298_STAGE4)
    im.io.write_npz(p.sub_problem(p.pedigree.neighbors(1298, depth=4, genotyped=False, add_parents=True)), im.itu.NBHRS1298_STAGE4_WITH_PARENTS)
    im.io.write_npz(p.sub_problem_of_parents(7, 462), im.itu.FAMILY13_STAGE4)
    im.io.write_npz(p.sub_problem_of_family(p.find_families_by_father(1715, genotyped=False)[0]), im.itu.FAMILY963_STAGE4)
    im.io.write_npz(p.sub_problem_of_parents(79, 324), im.itu.FAMILY_TOO_ZEROED_STAGE4)

    p = hutt('hutt.phased.npz')
    im.io.write_npz(p.sub_problem_of_parents(79, 324), im.itu.FAMILY_TOO_ZEROED_STAGE5)
    
def generate_test_sets_segments():
    '''Re-generate all chr22 IBD segment graph test sets.'''     
    # Generate IBD index. Requires phasing result + segment file for chr22 under the output directory
    im.index.index_segments.main(None, im.examples.CHR22 + '/hutt.phased.info.npz', im.examples.CHR22 + '/segments.out', im.itu.SEGMENT_INDEX_CHR22, num_processes=1, region_size=10, snp_index=1200, min_degree=4, min_length=4.0, margin=0.8, debug=2)

def single_validation(fraction=None, test_index=None):
    '''Run a single validation experiment with fraction% of deleted test genotypes. Returns a validation
    Experiment object.'''
    p = im.hutt('hutt.npz')
    e = im.v.Experiment(p, fraction=fraction, test_index=test_index)
    phaser = im.phase_main.main_phaser(print_times=True)
    e.run(phaser)
    (_, stats) = im.plots.plot_experiment_stats(e)
    print stats[np.argsort(-stats[:, 2] / stats[:, 1]), :]
    return e

def parametric_validation(input_data_set):
    '''Run the standard parametric validation on a NPZ data set under the OBER data directory.
    Returns a summary-statistics array.''' 
    p = im.io.read_npz(PREFIX + '/data/' + input_data_set)
    num_experiments = 10
    fraction = np.linspace(0.01, 0.05, num_experiments)
    phaser = im.phase_main.main_phaser()
    results = im.v.parametric_experiment(p, fraction, phaser)
    im.plots.plot_parametric_experiment(results)
    return results

#---------------------------------------------
# Methods - CGI IDs
#---------------------------------------------
def write_id_coefs_legacy(p, out_file, params=im.PhaseParam()):
    '''Given a legacy Problem object p, import ID coefficients from global ID coefficient file.'''   
    i = p.pedigree.genotyped_sample_id()
    a = np.array([(x, y, m[0]) + tuple(m[1].tolist()) for (x, y, m) in  ((x, y, params.id_coefs(x, y)) for x, y in itertools.product(i, i))])
    np.savetxt(out_file, a, fmt='%d %d %e %e %e %e %e %e %e %e %e %e')

'''Return the list of WGS Hutterities FINDIVs (sample ids).'''
wgs_sample_id = lambda sort_ids = True: db_gene.cgi.ids.sample_ids(sort_ids=sort_ids)

__wgs_sample_index = None
def wgs_sample_index():
    '''Return the list of WGS Hutterities pedigree sample indices.'''
    global __wgs_sample_index
    if __wgs_sample_index is None:
        ped = im.io_pedigree.read(im.itu.HUTT_PED, CHR22 + '/hutt.tfam')
        __wgs_sample_index = np.array(map(ped.node_of.get, db_gene.cgi.ids.sample_ids()))
    return __wgs_sample_index

__imputed_sample_id = None
def imputed_sample_id():
    '''Return the list of non-WGS Hutterities FINDIVs (sample ids).'''
    global __imputed_sample_id
    if __imputed_sample_id is None:
        ped = im.io_pedigree.read(im.itu.HUTT_PED, CHR22 + '/hutt.tfam')
        __imputed_sample_id = ped.sample_ids(im.imputed_sample_index())
    return __imputed_sample_id

__imputed_sample_index = None
def imputed_sample_index():
    '''Return the list of non-WGS Hutterities pedigree sample indices.'''
    global __imputed_sample_index
    if __imputed_sample_index is None:
        ped = im.io_pedigree.read(im.itu.HUTT_PED, CHR22 + '/hutt.tfam')
        __imputed_sample_index = np.array(list(set(xrange(ped.num_genotyped)) - set(wgs_sample_index())))
    return __imputed_sample_index

__sample_id_to_index = None
def sample_id_to_index():
    '''Return a sample FINDIV-to-index dictionary'''
    global __sample_id_to_index
    if __sample_id_to_index is None:
        ped = im.io_pedigree.read(im.itu.HUTT_PED, CHR22 + '/hutt.tfam')
        __sample_id_to_index = ped._node_of
    return __sample_id_to_index

__sample_index_to_id = None
def sample_index_to_id():
    '''Return a sample index-to-FINDIV-to dictionary'''
    global __sample_index_to_id
    if __sample_index_to_id is None:
        ped = im.io_pedigree.read(im.itu.HUTT_PED, CHR22 + '/hutt.tfam')
        __sample_index_to_id = ped._sample_id
    return __sample_index_to_id
