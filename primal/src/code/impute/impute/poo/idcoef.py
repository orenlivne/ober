#!/usr/bin/env python
'''
============================================================
Estimate detailed identity coefficients from IBD
cliques. Actually, returns counts of each identity state
for each pair of individuals. 

Created on March 12, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
# After flattening to n^2 x , t_counter corresponds to
# i[0] i[0] ... deltas ...
# i[0] i[1] ... deltas ...
# ... 
# where i=im.sample_index_to_id() or can be generally read from the corresponding phasing output.

import impute as im, os, numpy as np, time, sys, optparse, util, itertools as it

#---------------------------------------------
# Constants
#---------------------------------------------

# Edge bit constants representing the IBD state of two alleles (0,1) in two individuals (A,B)
# Each constants represents the existence (or lack of) of an edge in the 2x2 identity state graph.
A0_A1 = 1
A0_B0 = 2
A0_B1 = 4
A1_B0 = 8
A1_B1 = 16
B0_B1 = 32

# Detailed identity coefficients
S_star = np.array([
A0_A1 | A0_B0 | A0_B1 | A1_B0 | A1_B1 | B0_B1,
A0_A1 | A0_B0 | A1_B0,
A0_A1 | A0_B1 | A1_B1,
A0_B0 | A0_B1 | B0_B1,
A1_B1 | A1_B0 | B0_B1,
A0_A1 | B0_B1,
A0_A1,
B0_B1,
A0_B0 | A1_B1,
A0_B0,
A1_B1,
A0_B1 | A1_B0,
A0_B1,
A1_B0,
0])
T = np.zeros((64,), dtype=int)
T[S_star] = np.arange(len(S_star))

# Condensed identity state corresponding to each detailed identity state
S = np.array([0, 2, 2, 4, 4, 1, 3, 5, 6, 7, 7, 6, 7, 7, 8])

#---------------------------------------------
# Methods
#---------------------------------------------
def __parse_command_line_args(argv):
    '''Parse and validate command-line arguments.'''
    PROGRAM = os.path.basename(argv[0])
    usage = 'Usage: %s [flags] <segment-index> <out-file>\n\n' \
        'Return detailed identity state appearance counts from a segment\n' \
        'index of cliques.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-a', '--algorithm', type='str', dest='algorithm', default='full',
                      help='Imputation algorithm: sparse|full.')
    parser.add_option('-c', '--chr', type='int', dest='chrom', default=0,
                      help='Chromosome number to process. Overrides other options if non-zero. 0=do chromosomes in the range specified by the start, end chromosomes')
    parser.add_option('-r', '--region', type='int', dest='region', default= -1,
                      help='Only process a specific region index. If negative, all regions are processed')
    options, args = parser.parse_args(argv[1:]) 
    if len(args) != 2:
        print('\nMust specify two mandatory arguments. Found %d' % (len(args),))
        print usage
        sys.exit(1)
    if options.chrom < 1 or options.chrom > 22:
        print usage
        print('\nMust specify a chromosome number in 1..22.')
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
    return args, options

def num_regions(num_snps, region_size):
    '''Return the number of chromosomal regions given the total # of snps and a region size.'''
    return (num_snps / region_size) if num_snps % region_size == 0 else (num_snps / region_size + 1)  

def idcoef_sparse(ibd, chrom, regions): 
    '''Main call that calculates identity coefficients (if out='idcoef') or identity state counts
    (out='count') for chromosome chrom from the IBD segment index ibd. If region >= 0, only
    this IBD index region is processed, otherwise all regions are.'''  
    # Global metadata - SNPs
    ibd._load_chrom(chrom)
    bp = ibd._snp['base_pair']
    
    # Global metadata - samples
    n = ibd.num_samples
    samples = np.arange(n)

    # Counts #times each identity state between each two pairs of individuals occurs.
    # We only translate from 64 bits to 15 identity states at the end (faster but larger storage)
    t_counter = np.zeros((n, n, len(S_star)), dtype=np.uint32)
    # t = np.zeros((n, n), dtype=np.byte)  # Encodes identity state between each two pairs of individuals.
    total_bp = 0
    for region in regions:
        ibd._load_region_by_number(region)
        start_time = time.time()
        for snp, snp_cliques in enumerate(ibd._groups, ibd._start):  # [ibd._groups[99]]:  # ibd._groups
            # For each segment [bp[snp-1],bp[snp]], mark IBD state and increment t_count
            segment_len = 0 if (snp == 0) else (bp[snp] - bp[snp - 1])
            # t.fill(0)
            num_cliques = len(snp_cliques) - 1
            clique_comp = np.empty((num_cliques,), dtype=object)  # clique_comp[c] = list of samples NOT in clique #(c-1)
            hbd = np.empty((num_cliques,), dtype=object)  # hbd[c] = list of HBD samples in clique #(c-1)
            all_hbd = np.empty((0,), dtype=int)  # The union of hbd[c] for all c
            all_clique = np.empty((0,), dtype=int)  # The union of all clique samples
            for c, clique in enumerate(snp_cliques[1:]):
                print '\t', 'clique', (c + 1), 'size', len(clique)
                # Work with the smaller t-sub-matrix u of the samples in the 'sample' array
                sample, allele = clique[:, 0], clique[:, 1]
                clique_comp[c] = np.setdiff1d(samples, sample)
                all_clique = np.concatenate((all_clique, sample))
                m = len(sample)  # Size of clique
                u = np.zeros((m, m), dtype=np.byte)  # Encodes identity state between each two pairs of individuals.
                
                # Both (A,0),(B,0) are in clique ==> A0 = B0 IBD 
                i = np.where(allele == 0)[0]
                ind = np.meshgrid(i, i)
                u[ind] = u[ind] | A0_B0
     
                # Both (A,1),(B,1) are in clique ==> A1 = B1 IBD
                j = np.where(allele == 1)[0]
                ind = np.meshgrid(j, j)
                u[ind] = u[ind] | A1_B1
     
                # Both (A,0),(A,1) are in clique ==> mark A0 = A1 (row A) and B0 = B1 (column A)
                # This is only for samples within the clique; for HBD samples vs. samples outside
                # the clique, we update t-entries outside this loop
                ind = np.in1d(sample[i], sample[j])
                k = sample[i[ind]]
                hbd[c] = k
                all_hbd = np.concatenate((all_hbd, k))
                u[ind, :] = u[ind, :] | A0_A1
                u[:, ind] = u[:, ind] | B0_B1

                # (A,0),(B,1) are in clique and A != B ==> mark A0 = B1.
                # Similarly for the symmetric case (A,1), (B,0) ==> mark A1 = B0.
                ind = np.meshgrid(i, j)
                u[ind] = u[ind] | A0_B1
                ind = np.meshgrid(j, i)
                u[ind] = u[ind] | A1_B0
                
                # Increment counts of the identity states observed at this SNP for clique samples
                B, A = np.meshgrid(sample, sample)
                nz = (u != 0)
                t_counter[A[nz], B[nz], T[u[nz]]] += segment_len

            # Now, for all HBD samples vs. non-clique samples, mark A0 = A1 (row A) and B0 = B1 (column A)
            row_code, col_code, both_code = T[A0_A1], T[B0_B1], T[A0_A1 | B0_B1]
            
            for i, k in it.ifilter(lambda (_, k): k.size > 0, it.izip(clique_comp, hbd)):
                # HBD samples vs. non-clique non-HBD samples
                I = np.setdiff1d(i, all_hbd)
                p, q = np.meshgrid(k, I)
                t_counter[p, q, row_code] += segment_len
                p, q = np.meshgrid(I, k)
                t_counter[p, q, col_code] += segment_len
                # HBD samples vs. non-clique HBD samples
                I = np.intersect1d(i, all_hbd)
                p, q = np.meshgrid(k, I)
                t_counter[p, q, both_code] += segment_len
                p, q = np.meshgrid(I, k)
                t_counter[p, q, both_code] += segment_len
            
            # Haplotype vs. itself (not included in the logic above this line to avoid double-counting)
            t_counter[samples, samples, T[A0_B0 | A1_B1]] += segment_len
            
            total_bp += segment_len
            print region, snp, bp[snp], segment_len, total_bp 
        print 'Region', region, '#cliques', len(snp_cliques[1:]), 'time', time.time() - start_time
    return t_counter, total_bp

def append_to_nz_list(nz, x):
    return tuple(np.concatenate((nz[i], x[i].flatten()), axis=0) for i in xrange(len(nz)))

def idcoef_sparse2(ibd, chrom, regions): 
    '''Main call that calculates identity coefficients (if out='idcoef') or identity state counts
    (out='count') for chromosome chrom from the IBD segment index ibd. If region >= 0, only
    this IBD index region is processed, otherwise all regions are.
    
    Slow implementation that updates the full t_counter matrix and requires x3 more storage.
    For debugging purposes only. Should yield identical results to idcoef_sparse.
    
    Slower than idcoef_full()...'''  
    # Global metadata - SNPs
    ibd._load_chrom(chrom)
    bp = ibd._snp['base_pair']
    
    # Global metadata - samples
    n = ibd.num_samples
    samples = np.arange(n)
    B, A = np.meshgrid(np.arange(n), np.arange(n))

    # Counts #times each identity state between each two pairs of individuals occurs.
    # We only translate from 64 bits to 15 identity states at the end (faster but larger storage)
    t_counter = np.zeros((n, n, len(S_star)), dtype=np.uint32)
    t = np.zeros((n, n), dtype=np.byte)  # Encodes identity state between each two pairs of individuals.
    total_bp = 0
    for region in regions:
        ibd._load_region_by_number(region)
        start_time = time.time()
        for snp, snp_cliques in enumerate(ibd._groups, ibd._start):  # [ibd._groups[99]]:  # ibd._groups
            # For each segment [bp[snp-1],bp[snp]], mark IBD state and increment t_count
            segment_len = 0 if (snp == 0) else (bp[snp] - bp[snp - 1])

            # Segments of a haplotype an itself (not included in our cliques/segment files)
            t.fill(0)
            t[samples, samples] = t[samples, samples] | (A0_B0 | A1_B1)
            nz = (samples, samples)  # Keeps track of the indices of all non-zeros in t
            
            for _, clique in enumerate(snp_cliques[1:]):
                # print '\t', 'clique', (c+1), 'size', len(clique) 
                sample, allele = clique[:, 0], clique[:, 1]
                
                # Both (A,0),(B,0) are in clique ==> A0 = B0 IBD 
                i = sample[allele == 0]
                ind = np.meshgrid(i, i)
                t[ind] = t[ind] | A0_B0
                nz = append_to_nz_list(nz, ind)
    
                # Both (A,1),(B,1) are in clique ==> A1 = B1 IBD
                j = sample[allele == 1]
                ind = np.meshgrid(j, j)
                t[ind] = t[ind] | A1_B1
                nz = append_to_nz_list(nz, ind)
    
                # Both (A,0),(A,1) are in clique ==> HBD; mark A0 = A1 (row A) and B0 = B1 (column A)
                k = np.intersect1d(i, j)
                t[k, :] = t[k, :] | A0_A1 
                t[:, k] = t[:, k] | B0_B1
                nz = append_to_nz_list(nz, np.meshgrid(k, samples))
                nz = append_to_nz_list(nz, np.meshgrid(samples, k))
    
                # (A,0),(B,1) are in clique and A != B ==> mark A0 = B1.
                # Similarly for the symmetric case (A,1), (B,0) ==> mark A1 = B0.
                ind = np.meshgrid(i, j)
                t[ind] = t[ind] | A0_B1
                nz = append_to_nz_list(nz, ind)

                ind = np.meshgrid(j, i)
                t[ind] = t[ind] | A1_B0
                nz = append_to_nz_list(nz, ind)
                
            # Increment identity states counts with those observed at this SNP
            # Only increment non-zero values
            t_counter[A[nz], B[nz], T[t[nz]]] += segment_len
            total_bp += segment_len
            print region, snp, bp[snp], segment_len, total_bp
            
        print 'Region', region, '#cliques', len(snp_cliques[1:]), 'time', time.time() - start_time

    return t_counter, total_bp

def idcoef_full(ibd, chrom, regions): 
    '''Main call that calculates identity coefficients (if out='idcoef') or identity state counts
    (out='count') for chromosome chrom from the IBD segment index ibd. If region >= 0, only
    this IBD index region is processed, otherwise all regions are.
    
    Slow implementation that updates the full t_counter matrix and requires x3 more storage.
    For debugging purposes only. Should yield identical results to idcoef_sparse.'''  
    # Global metadata - SNPs
    ibd._load_chrom(chrom)
    bp = ibd._snp['base_pair']
    
    # Global metadata - samples
    n = ibd.num_samples
    samples = np.arange(n)
    B, A = np.meshgrid(np.arange(n), np.arange(n))
    B, A = B[np.newaxis], A[np.newaxis]

    # Counts #times each identity state between each two pairs of individuals occurs.
    # We only translate from 64 bits to 15 identity states at the end (faster but larger storage)
    t_counter = np.zeros((n, n, 64), dtype=np.uint32)
    t = np.zeros((n, n), dtype=np.byte)  # Encodes identity state between each two pairs of individuals.
    total_bp = 0
    for region in regions:
        ibd._load_region_by_number(region)
        start_time = time.time()
        for snp, snp_cliques in enumerate(ibd._groups, ibd._start):  # [ibd._groups[99]]:  # ibd._groups
            # For each segment [bp[snp-1],bp[snp]], mark IBD state and increment t_count
            segment_len = 0 if (snp == 0) else (bp[snp] - bp[snp - 1])

            t.fill(0)
            # Segments of a haplotype an itself (not included in our cliques/segment files)
            t[samples, samples] = t[samples, samples] | (A0_B0 | A1_B1)
            for _, clique in enumerate(snp_cliques[1:], 1):
                # print '\t', 'clique', c, 'size', len(clique) 
                sample, allele = clique[:, 0], clique[:, 1]
                
                # Both (A,0),(B,0) are in clique ==> A0 = B0 IBD 
                i = sample[allele == 0]
                ind = np.meshgrid(i, i)
                t[ind] = t[ind] | A0_B0
    
                # Both (A,1),(B,1) are in clique ==> A1 = B1 IBD
                j = sample[allele == 1]
                ind = np.meshgrid(j, j)
                t[ind] = t[ind] | A1_B1
    
                # Both (A,0),(A,1) are in clique ==> HBD; mark A0 = A1 (row A) and B0 = B1 (column A)
                k = np.intersect1d(i, j)
                t[k, :] = t[k, :] | A0_A1 
                t[:, k] = t[:, k] | B0_B1 
    
                # (A,0),(B,1) are in clique and A != B ==> mark A0 = B1.
                # Similarly for the symmetric case (A,1), (B,0) ==> mark A1 = B0.
                ind = np.meshgrid(i, j)
                t[ind] = t[ind] | A0_B1
                ind = np.meshgrid(j, i)
                t[ind] = t[ind] | A1_B0
                
            # Increment counts of the identity states observed at this SNP
            t_counter[A, B, t] += segment_len
            total_bp += segment_len
            print region, snp, bp[snp], segment_len, total_bp
        print 'Region', region, '#cliques', len(snp_cliques[1:]), 'time', time.time() - start_time

    t_counter = t_counter[:, :, S_star]  # .astype(np.double) / snps
    return t_counter, total_bp

def idcoef(ibd, chrom, algorithm='sparse', out='idcoef', do_region= -1):
    if algorithm == 'sparse': f = idcoef_sparse2
    elif algorithm == 'full': f = idcoef_full
    else: raise ValueError('Unsupported algorithm ''%s''' % (algorithm,))
    
    regions = [do_region] if do_region >= 0 else xrange(num_regions(len(ibd._snp), ibd._region_size))
    t_counter, total_bp = f(ibd, chrom, regions)
     
    # Save counts
    if out == 'idcoef':
        return t_counter.astype(np.double) / total_bp
    else:
        return t_counter, np.int64(total_bp)

def reshape_idcoef_array(t_counter):
    '''Reshape a 3-D idcoef array to a 2-D.'''
    n = t_counter.shape[0]
    return t_counter.reshape((n * n, t_counter.shape[2]))

####################################################################################
if __name__ == '__main__':
    args, options = __parse_command_line_args(sys.argv)
    ibd = im.index.segment_index.SegmentIndex(args[0])
    t_counter, total_bp = idcoef(ibd, options.chrom, algorithm=options.algorithm, out='count', do_region=options.region)
    # Save counts
    np.savetxt(args[1], reshape_idcoef_array(t_counter), fmt='%d')
