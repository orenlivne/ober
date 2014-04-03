#!/usr/bin/env python
'''
============================================================
Run manual imputation on monogenic SNPs (Jessica's rare
SNPs) for which there were discordances between the CGI and
Daokta genotypes. 

Created on June 12, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os, util, sys, numpy as np
from impute.ibd.index import index_segments_beagle

def pipeline_monogenic_validation(work_dir=os.environ['OBER_OUT'] + '/requests/monogenic/work',
                                  index_segments_dir=os.environ['OBER_OUT'] + '/requests/monogenic/work/index_segments',
                                  region_size=100,
                                  theta_affinity=0.95,
                                  theta_weight=0.5,
                                  regenerate_segments=True,
                                  snps=None,  # np.array([6, 8]),
                                  debug=1,
                                  debug_sample=512):
    # Load SNPs
    problem = im.io.read_plink(prefix=work_dir + '/monogenic.12', pedigree=im.itu.HUTT_PED, haplotype=None, frames=None)
    # Testing: simulate aligned samples output (hap types should be 2 in the imputed genotype output line)
    problem.haplotype.poo_phase = np.zeros((problem.num_samples,), dtype=np.byte)
    problem.haplotype.poo_phase[np.array([0, 1])] = 1
    problem.haplotype.poo_phase[np.array([2, 3])] = -1
    
    # Create segments only for the regions around each snp
    if regenerate_segments:
        for row in (problem.info.snp[snps] if snps is not None else problem.info.snp):
            # Find SNP's region (the one containing its base-pair position) 
            chrom, bp = row['chrom'], row['base_pair']
            phasing_dir = '%s/phasing/chr%d' % (os.environ['OBER_OUT'], chrom)
            index_segments_chrom_dir = '%s/chr%d' % (index_segments_dir, chrom)
            info_file = '%s/hutt.phased.info.npz' % (phasing_dir,)
            info = im.io.read_info_npz(info_file)
            snp_bp = info.snp['base_pair']
            snp_index = util.nearest_neighbor_in_list_tree(bp, snp_bp, util.list_index_tree(snp_bp))
            snp_index = snp_index if snp_bp[snp_index] <= bp else snp_index - 1
            start = region_size * (snp_index / region_size)
            stop = start + region_size
            segment_file = '%s/segments-%d-%d.out' % (index_segments_chrom_dir, start, stop)
            if not os.path.exists(segment_file):
                util.mkdir_if_not_exists(index_segments_chrom_dir)
                util.run_command('find-segments-of-snp-range %d %d < %s/segments.out > %s' % (start, stop, phasing_dir, segment_file)) 
            
            # Index segments
            if regenerate_segments or \
            not os.path.exists('%s/metadata.npz' % (index_segments_chrom_dir,)) or \
            not os.path.exists('%s/region-%d.npz' % (index_segments_chrom_dir, start)):
                index_segments_beagle.main(segment_file, info_file, segment_file, index_segments_chrom_dir,
                                           snp_index=snp_index, debug=2,
                                           theta_affinity=theta_affinity, theta_weight=theta_weight)
    
    # Impute using the newly generated segment index
    _, t = im.v.iv.impute_problem(problem, debug=debug, remove_partial_calls=True,
                                  segment_location=index_segments_dir,  # if regenerate_segments else None,
                                  snps=snps, debug_sample=debug_sample)

    im.io.write_plink(im.Problem(genotype=t.imputed, pedigree=im.examples.hutt_pedigree(), haplotype=None, frames=None),
                      work_dir + '/imputed.12', save_frames=False, save_haplotype=False)
    im.cgi.io_cgi.write_imputed(t, sys.stdout, poo_phase=problem.haplotype.poo_phase)
    with open(work_dir + '/imputed.12.lgen', 'wb') as f:
        im.cgi.io_cgi.write_imputed_lgen(t, f)
    return t
