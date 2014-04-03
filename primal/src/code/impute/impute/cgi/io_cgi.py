'''
============================================================
Read and write whole-genome sequencing CGI data sets
and our imputed data sets in a similar format.

Created on July 25, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, StringIO, itertools as it
from impute.tools.recode import CGI_LETTER_TO_ALLELE
from util import flattened_meshgrid

#---------------------------------------------
# Methods
#---------------------------------------------
def read_genotype(chrom, sample_id, index, lines):
    '''Read CGI genotypes from a list of lines.'''
    #-------------------------------------------
    # Allocate genotype data arrays
    #-------------------------------------------
    num_snps = len(lines)
    num_samples = len(sample_id)
    # Standard SNP metadata
    snp = np.zeros((num_snps,),
                   dtype=[('chrom', np.uint8),  # Chromosome # containing the SNP
                          ('name', np.chararray),  # SNP name (e.g., 'rs...')
                          ('dist_cm', np.float),  # Genetic distance from beginning of chomorsome
                          ('base_pair', np.uint)  # Base pair position on chromosome
                          ])
    # All SNP attributes (allele letters, RS #, etc.)
    metadata = [[None]] * num_snps
    # Sample genotypes 
    data = np.zeros((num_snps, num_samples, 2), dtype=np.byte)

    #-------------------------------------------
    # Bulk-extract data of the selected lines
    #-------------------------------------------
    for row, line in enumerate(lines):
        snp[row] = (chrom, line[7], 0, int(line[3]))
        metadata[row] = line[0:8]
        data[row, :] = np.array(sum(([CGI_LETTER_TO_ALLELE[item[0]], CGI_LETTER_TO_ALLELE[item[1]]]
                                     for item in line[8:]), [])).reshape(num_samples, 2)[index, :]
#        d = data[row, :]
#        _, f2 = im.gt.allele_frequencies(d)
#        if f2 > 0.7:  # Minor allele has a much larger frequency than the major allele, swap them
#            #print 'SNP chr%d:%d-%d (%s): Major-minor allele swap detected, fixing' % \
#            #(chrom, int(line[2]), int(line[3]), line[7])
#            im.gt.swap_alleles(d)
    g = im.factory.GenotypeFactory.new_instance('genotype', data, snp, sample_id=sample_id)
    g.metadata = metadata
    return g

def write_imputed(t, out, debug=False, poo_phase=None):
    '''Write imputed genotypes to the stream out in CGI format.'''
    data, metadata, hap_type = im.recode.recode_cgi(t.imputed_data), t.genotype.metadata, t.imputed_hap_type
    if poo_phase is not None:
        aligned_samples = np.where(poo_phase)[0]
        t.imputed_hap_type[:, aligned_samples] = im.constants.PHASED_WITH_ORIGIN
        # Flip haplotypes of samples with flipped POO phase
        flipped_samples = np.where(poo_phase < 0)[0]
        orig = flattened_meshgrid(flipped_samples, im.constants.ALLELES)
        flipped = flattened_meshgrid(flipped_samples, list(reversed(im.constants.ALLELES)))
        data[:, orig[0], orig[1]] = data[:, flipped[0], flipped[1]] 
        
    if debug: np.set_printoptions(threshold=np.nan)
    for snp in t.genotype.snp_range:
        # Ensure that all fields are non-empty - easier to parse by subsequent processes
        if metadata:
            np.savetxt(out, np.array(map(lambda x: x if x else '-', metadata[snp])), fmt='%s', newline='\011', delimiter='')
        # Remove trailing tab at the end of the line produced by the numpy savetxt call
        out_str = StringIO.StringIO()
        np.savetxt(out_str, [x[0] + x[1] for x in it.izip(it.imap(str, hap_type[snp]), (g[0] + g[1] for g in data[snp]))], fmt='%s', newline='\011', delimiter='')
        out.write(out_str.getvalue()[:-1])
        out.write('\n')
        out.flush()
        if debug: print np.concatenate((np.arange(data.shape[1])[np.newaxis].transpose(), hap_type[snp][np.newaxis].transpose(), t.imputed_data[snp]), axis=1) 
    if debug: np.set_printoptions(threshold=1000)

def write_imputed_lgen(t, out, debug=False):
    '''Write imputed genotypes to the stream out in PLINK LGEN format.'''
    data, snp_names, sample_id = im.recode.recode_cgi(t.imputed_data), t.genotype.snp['name'], t.pedigree.sample_id
    for snp in t.genotype.snp_range: out.write('\n'.join('%s\t%d\t%s\t%s' % (snp_names[snp], sample_id[i], g[0], g[1]) for i, g in enumerate(data[snp])) + '\n')

#---------------------------------------------
# Private Methods
#---------------------------------------------
