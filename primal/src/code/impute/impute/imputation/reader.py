#!/usr/bin/env python
'''
============================================================
Functions and classes that load data from related files.

Created on December 3, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, csv, glob, os, db_gene
from impute.ibd.segment import START, STOP
from impute.data.factory import GenotypeFactory
from db_gene.cgi.ids import DEFAULT_ID_INDEX

#---------------------------------------------
# Methods - Affymetrix phasing Problem reading
#---------------------------------------------
def problem_to_imputation_set(p, index_file=DEFAULT_ID_INDEX, genotypes_as_is=False):
    '''Create a training set by extracting the training sample affymetrix genotypes from a Problem object.'''
    training_id = db_gene.cgi.ids.cgi_id_to_sample_id(file_name=index_file).values()  # @UndefinedVariable
    if genotypes_as_is: 
        g = im.factory.GenotypeFactory.new_instance('genotype', p.g.copy(), p.genotype.snp, sample_id=training_id)
    else:
        training_index = np.array([p.pedigree.node_of[x] for x in training_id])
        g = im.factory.GenotypeFactory.new_instance('genotype', p.g[:, training_index, :].copy(), p.genotype.snp, sample_id=training_id)
    return im.imputation.ImputationSet(p.pedigree, g)

#---------------------------------------------
# Methods - iPlex Genotype Reading
#---------------------------------------------
def iplex_to_genotype(path, t, suffix='_iPlex_HUTTERITE.txt'):
    '''Convert an iPlex data into a Genotype object that matches the ordering of SNPs and
    samples in an ImputationSet object t.'''
    # Read iPlex data into a dictionary of dictionaries
    iplex = __read_iplex(path, suffix)
    
    # Match SNPs by name. This is a slow (m x n) operation where m=#iplex SNPs, n=#imputation_set SNPs,
    # but speed is not an issue here since m, n are in the tens to hundreds
    iplex_snps = iplex.keys()
    print iplex_snps
    snp_metadata = t.snp['name']
    # index of each IPLEX SNP in imputation_set
    index = [np.where([iplex_snp in x for x in snp_metadata])[0][0] for iplex_snp in iplex_snps]
    
    # Match genotype IDs of each iplex SNP and imputation_set
    s = t.pedigree.genotyped_sample_id()
    data = np.zeros_like(t.imputed_data)
    print t.genotype.map
    num_snps = len(iplex_snps)
    print 'num_snps', num_snps
    for i in xrange(num_snps):
        snp = iplex_snps[i]
        snp_index = index[i]
        iplex_data = iplex[snp]
        values = np.array(iplex_data.values())
        
        letter = t.genotype.map[snp_index]
        print '#%2d iPlex SNP %-14s index in t %3d %s/%s' % (i, snp, snp_index, letter[0], letter[1])
        # Reverse strand if needed
        if values[0][0] not in letter:
            print 'Reversing strand letters of iPlex data (Dakota labeling uses standard dbSNP letters)'
            letter = ''.join(REVERSE_STRAND[x] for x in letter)

        allele = __allele_dict(letter)
        iplex_id = np.array(iplex_data.keys())
        samples_in_t = np.in1d(iplex_id, s)
        sample_index = [t.pedigree.node_of[x] for x in iplex_id[samples_in_t]]

        recoded_genotypes = np.array([[allele[v[0]], allele[v[1]]] for v in values[samples_in_t]])
        f1, f2 = im.gt.allele_frequencies(recoded_genotypes)
        if f2 > 0.7:  # Minor allele has a much larger frequency than the major allele, swap them
            print 'Major-minor allele CGI letter swap detected, fixing (f1=%.2f, f2=%.2f)' % (f1, f2)
            im.gt.swap_alleles(recoded_genotypes)
            t.genotype.map[snp_index] = letter[1] + letter[0]
        
        # Insert data into the appropriate row of the target genotype object data array 
        data[snp_index, sample_index, :] = recoded_genotypes
    # return data, snp_metadata[index], t.genotype.sample_id
    g = GenotypeFactory.new_instance('genotype', data, t.snp, t.genotype.sample_id)
    g.map = t.genotype.map
    return g 

#---------------------------------------------
# Private Methods
#---------------------------------------------        
'''Base-pair letter on the reverse strand.'''
REVERSE_STRAND = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

'''Read a single SNP''s iPlex data from a file into a dictionary of ID-to-genotypes.
Skip first line; extract non-control samples only.'''
__read_iplex_file = lambda f: dict((int(line[0]), line[1] + line[2]) for i, line in enumerate(csv.reader(f, delimiter='\t', skipinitialspace=True)) if i >= 1 and line and len(line) == 3 and line[0] != 'CONTROL')

'''Convert the two allele letters to a dictionary that maps iPlex symbols to allele numerical coding.'''
__allele_dict = lambda letter: {'0': 0, letter[0]: 1, letter[1]: 2}

'''Read iPlex genotypes (genotyped by Kristen; extracted from
the dakota database by William) from files in the current
directory, and convert to a dictionary of dictionaries, one per SNP. Each SNP's value is
a dictionary of ID-to-genotypes, as different SNPs may have different samples.

path = directory containing all iplex txt files
suffix = file suffix to look for under path'''
__read_iplex = lambda path, suffix: dict(((os.path.basename(f)[:-len(suffix)], __read_iplex_file(open(f, 'rb')))  for f in glob.glob(path + '/*' + suffix)))
        
####################################################################################
class PhasingResultReader(object):
    '''Reads metadata and retrieves haplotypes from a phasing output data directory.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, phasing_dir, name):
        '''Initialize a reader to read the data set named ''name'' (=filename prefix)
        from the phasing output directory ''phasing_dir''.'''
        # One underscore = protected member
        self.__phasing_dir = phasing_dir        
        self.__name = name
        phasing_stats_file_name = '%s/reduce/result/%s.stats.npz' % (phasing_dir, name)
        self.pedigree = np.load(phasing_stats_file_name)['pedigree'][0]

        chr_endpoints_file_name = phasing_dir + '/part_count.range.npz'
        a = np.load(chr_endpoints_file_name)
        self.__part_size = a['part_size'][0]
        self.__chrom_endpoints = a['endpoints']
        self.__part_count = np.array([int(np.ceil((x[STOP] - x[START]) / self.__part_size)) for x in self.__chrom_endpoints])
        # print self.__chrom_endpoints
        # print self.__part_count
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------        
    def get_window(self, chrom, bp, min_ibs_len_snp):
        '''Returns a haplotype matrix n x m x 2, where m = # samples and n is the number of SNPs,
        which include at most min_ibs_len_snp SNPs to the left and min_ibs_len_snp SNPs to right of the
        base pair position on chromosome chrom. Assumes that min_ibs_len_snp <= chrom part size. Under
        the hood, may retrieve from one or two part files and patch them together.
        
        Also returns the relative position of the nearest-right-neighbor of the bp within the window.'''
        part = self.__part_num(chrom, bp)
        this_part = self.__load_part(chrom, part)
        # k = right-neighboring-SNP of bp
        k = this_part.genotype.nearest_snp(bp)
        if this_part.genotype.snp['base_pair'][k] < bp:
            k += 1
        left = k - min_ibs_len_snp + 1  # inclusive
        right = k + min_ibs_len_snp  # non-inclusive
        print 'bp', bp, 'absolute k', k
        if left < 0 and part > 0:
            print 'straddling prev,this'            
            # SNP straddled between prev_part, this_part (in particular, this_part is not the first part)
            prev_part = self.__load_part(chrom, part - 1)
            print part - 1, left + prev_part.num_snps, prev_part.num_snps
            print part, 0, right
            data = np.concatenate((prev_part.h[left + prev_part.num_snps:prev_part.num_snps, :, :],
                                   this_part.h[0:right, :, :]),
                                  axis=0)
            snps = np.concatenate((prev_part.genotype.snp[left + prev_part.num_snps:prev_part.num_snps],
                                   this_part.genotype.snp[0:right]))
        elif right >= this_part.num_snps and part < self.__part_count[chrom - 1] - 1:
            # SNP straddled between this_part, next_part (in particular, this_part is not the last part)
            print 'straddling this,next'
            next_part = self.__load_part(chrom, part + 1)
            print part, left, this_part.num_snps
            print part + 1, 0, right - this_part.num_snps
            data = np.concatenate((this_part.h[left:this_part.num_snps, :, :],
                                   next_part.h[0:right - this_part.num_snps, :, :]),
                                  axis=0)
            snps = np.concatenate((this_part.genotype.snp[left:this_part.num_snps],
                                   next_part.genotype.snp[0:right - this_part.num_snps]))
        else:
            # Assuming min_ibs_len_snp <= part_size, there must be no straddling in this case
            print 'No straddling, this_part.num_snps', this_part.num_snps
            left = max(left, 0)
            right = min(right, this_part.num_snps)
            print part, left, right
            data = this_part.h[left:right, :, :]
            snps = this_part.genotype.snp[left:right]
        return im.factory.GenotypeFactory.new_instance('haplotype', data, snps,
                                                       sample_id=self.pedigree.sample_id), k - left
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------        
    def __part_num(self, chrom, bp):
        '''Return the part number that a base-pair position lies in. Raises an exception if we are
        outside the range.'''
        (start, stop) = self.__chrom_endpoints[chrom - 1]
        if bp < start or bp > stop:
            raise ValueError('Base-pair position chr%d:%d is outside the covered range [%d,%d]' % \
                             (chrom, bp, start, stop))
        return int(np.floor((bp - 1 - self.__chrom_endpoints[chrom - 1][START]) / self.__part_size))

    def __load_part(self, chrom, part):
        '''Load a chromsome part file from the output directory.'''
        return im.io.read_npz('%s/map/chr%d/%s_chr%d_part%d.npz' % \
                              (self.__phasing_dir, chrom, self.__name, chrom, part)) 
