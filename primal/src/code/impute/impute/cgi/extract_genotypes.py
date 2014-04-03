#!/usr/bin/env python
'''
============================================================
Convert Hutterites WGS data at SNP locations read from
standard input to PLINK TPED format.

Prerequisites: tabix

Created on November 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, csv, util, db_gene, numpy as np, StringIO, gzip, impute as im
from impute.tools.recode import CGI_LETTER_TO_ALLELE
from impute.data.constants import CHROMOSOMES

#---------------------------------------------
# Constants
#---------------------------------------------
# Path to tabix program
_TABIX_CMD = 'tabix' 

# Convert chromosome numeric value to letter
_to_letter = dict([(i, str(i)) for i in CHROMOSOMES] + [(23, 'X'), (24, 'Y'), (25, 'M')])

#---------------------------------------------
# Methods
#---------------------------------------------
def extract_genotypes(input_file,
                      var_file_prefix=os.environ['OBER_DATA'] + '/cgi/all.2012-09-20.testvar.chr',
                      index_file=os.environ['OBER_DATA'] + '/cgi/README.assembly_sample_subject.csv',
                      allele_swap=False, debug=False, custom_id=False):
    '''Extract SNP data from tabixed CGI var files and convert to a Genotype object.
    Base-pair ranges are read from input_file in the format: chr start_bp end_bp
    * Multiple ranges can be specified (one per line).
    * If input_file = ''-'', input is taken from standard input.
    * Use 23 for chrX, 24 for chrY. (BUT: chrs X, Y, M are not yet supported. TBA.)
    * By default, allele labels are swapped so that 1 is always the major and 2 is always the minor.'''
    
    # Read CGI-ID-to-our-sample-ID (FINDIV) dictionary
    cgi_id_to_sample_id = db_gene.cgi.ids.cgi_id_to_sample_id(index_file)
    
    # Read SNPs from stdin; assuming numeric chromosome numbers (typically we only use autosomal)
    input_snps = [(int(line[0]), (int(line[1]), int(line[2])))
            for line in csv.reader(input_file, skipinitialspace=True, delimiter=' ')]
    num_snps = len(input_snps)
    
    # Holds the sorted list sample IDs; assumed to be the same for all chromosomes;
    # if not, raise an exception in the chromosome loop below
    sample_id = None
    
    # Holds SNP metadata
    snp = np.zeros((num_snps,),
                   dtype=[('chrom', np.uint8),  # Chromosome # containing the SNP
                          ('name', np.chararray),  # SNP name (e.g., 'rs...')
                          ('dist_cm', np.float),  # Genetic distance from beginning of chomorsome
                          ('base_pair', np.uint)  # Base pair position on chromosome
                          ])
    # Genetic map (allele letters of each SNP)
    genetic_map = np.zeros((num_snps,), dtype=np.object)
    
    # Sample genotypes at all SNPs 
    data = None

    # Bulk-extract data using tabix for each chromosome
    row = 0
    for chrom, bps in util.to_set_dict(input_snps).iteritems():
        chrom_letter = _to_letter[chrom]
        var_file = '%s%s.tsv.gz' % (var_file_prefix, chrom_letter)

        # Read header line from the compressed var file of this chromosome to get CGI IDs
        for line in csv.reader(gzip.GzipFile(var_file), skipinitialspace=True, delimiter='\t'):
            raw_id = [cgi_id_to_sample_id[cgi_id] for cgi_id in line[8:]]
            chr_sample_id, index = util.sort_with_index(np.array(raw_id))
            break
        if sample_id is None:
            # Save global IDs 
            sample_id = chr_sample_id
            num_samples = len(sample_id)
        elif not np.array_equal(chr_sample_id, sample_id):
            # Make sure the IDs are the same across all chromosomes
            raise ValueError('Sample IDs are not the same across all chromosomes:\n' 
                             + sample_id + '\n' + chr_sample_id)
        
        # Read SNP data from the compressed chromosome archive using tabix; append to genotype array
        if data is None:
            data = np.zeros((num_snps, num_samples, 2), dtype=np.byte)
        _, output, _ = util.run_command(_TABIX_CMD + ' ' + var_file + ' '
                                        + ' '.join('chr%s:%d-%d' % (chrom_letter, x[0], x[1]) for x in bps),
                                        verbose=debug)
        # Apparently there's a tabix bug: some records might not be contained in the requested region.
        # Only process lines that are within the requested regions
        relevant_lines = [a for a in csv.reader(StringIO.StringIO(output), skipinitialspace=True, delimiter='\t')
                          if np.array([(x[0] <= int(a[2])) & (int(a[3]) <= x[1]) and a[4] == 'snp'
                                       for x in bps]).any()]
        for line in relevant_lines:
            if len(snp) <= row:
                snp.resize(2 * len(snp))
                genetic_map.resize(2 * len(snp))
                data.resize((2 * len(snp), num_samples, 2), refcheck=False)
            # If SNP name is missing, create a unique identifier by concatenating the chromosome, begin and minor allele
            snp[row] = (chrom, line[7] if line[7] and not custom_id else 'chr%s_%d_%s' % (chrom_letter, int(line[2]), line[6]), 0, int(line[3]))
            genetic_map[row] = line[5] + line[6]
            # print line[8:]
            if chrom <= 22:
                # Autosomes - genotype data (two letters)
                data[row, :] = np.array(sum(([CGI_LETTER_TO_ALLELE[item[0]], CGI_LETTER_TO_ALLELE[item[1]]]
                                             for item in line[8:]), [])).reshape(num_samples, 2)[index, :]
            else:
                # X/Y chromosomes: if one letter is available => duplicate it to make the individual homozygous
                # (or 00 for a single missing value)
                data[row, :] = np.array(sum(([CGI_LETTER_TO_ALLELE[item[0]], CGI_LETTER_TO_ALLELE[item[1 if len(item) == 2 else 0]]]
                                             for item in line[8:]), [])).reshape(num_samples, 2)[index, :]
            d = data[row, :]
            _, f2 = im.gt.allele_frequencies(d)
            if allele_swap and f2 > 0.7:  # Minor allele has a much larger frequency than the major allele, swap them
                print 'chr%d:%d-%d (%s): Major-minor allele swap detected, fixing' % \
                (chrom, int(line[2]), int(line[3]), line[7])
                im.gt.swap_alleles(d)
            row += 1

    # Not all SNPs might be found; restrict data arrays to the first row rows 
    snp = snp[0:row]
    genetic_map = genetic_map[0:row].tolist()
    data = data[0:row]
    print 'Found %d SNP(s) in %d requested base-pair range(s)' % (row, num_snps)
    
    # Construct Genotype object; write in desired output format
    g = im.factory.GenotypeFactory.new_instance('genotype', data, snp, sample_id=sample_id)
    g.map = genetic_map
    return g
